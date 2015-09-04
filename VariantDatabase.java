package nhs.genetics.cardiff;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.neo4j.graphdb.*;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;

import java.io.File;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by ml on 23/06/15.
 */

//TODO convert VCF to minimal represented and remove sample-specific info; annotate then use as input
//TODO import dbsnfp and use e v79
//TODO parse OMIM records correctly
//TODO retrieve pubmed abstracts (web ui)
//TODO check alamut for extra functionality
//TODO skip NTC
//TODO add variant type label
//TODO add sample type

public class VariantDatabase {

    private static final Logger log = Logger.getLogger(VariantDatabase.class.getName());

    private File neo4jDBPath;
    private GraphDatabaseService graphDb;
    private VCFFileReader variantVcfFileReader, annotationVcfFileReader;
    private HashMap<String, HashSet<String>> geneMap2;
    private ArrayList<VariantContext> vcfBody = new ArrayList<>(); //VCF file body
    private HashMap<GenomeVariant, HashSet<VEPAnnotation>> vepAnnotations = new HashMap<>(); //all VEP annotations
    private HashMap<GenomeVariant, HashMap<String, Double>> populationFrequencies = new HashMap<>(); //population frequencies from mulitple sources
    private HashMap<GenomeVariant, Node> newVariantNodes = new HashMap<>(); //new variants added during this session
    private HashMap<GenomeVariant, Node> existingNodes = new HashMap<>();
    private HashMap<String, Node> runInfoNodes = new HashMap<>(); //samples added during this session
    private HashMap<GenomeVariant, HashMap<String, Node>> annotationNodes = new HashMap<>(); //
    private HashMap<String, Node> symbolNodes = new HashMap<>(); //symbols added during this session
    private HashMap<String, Node> featureNodes = new HashMap<>(); //features added during this session

    private Label sampleLabel = DynamicLabel.label("Sample");
    private Label variantLabel = DynamicLabel.label("Variant");
    private Label annotationLabel = DynamicLabel.label("Annotation");
    private Label symbolLabel = DynamicLabel.label("Symbol");
    private Label canonicalLabel = DynamicLabel.label("Canonical");
    private Label featureLabel = DynamicLabel.label("Feature");
    private Label disorderLabel = DynamicLabel.label("Disorder");
    private Label runInfoLabel = DynamicLabel.label("RunInfo");

    public VariantDatabase(VCFFileReader variantVcfFileReader, VCFFileReader annotationVcfFileReader, File neo4jDBPath, HashMap<String, HashSet<String>> geneMap2){
        this.variantVcfFileReader = variantVcfFileReader;
        this.neo4jDBPath = neo4jDBPath;
        this.annotationVcfFileReader = annotationVcfFileReader;
        this.geneMap2 = geneMap2;
    }

    private enum relTypes implements RelationshipType
    {
        HAS_HET_VARIANT,
        HAS_HOM_VARIANT,
        IN_SYMBOL,
        IN_FEATURE,
        HAS_UNKNOWN_CONSEQUENCE,
        HAS_ASSOCIATED_DISORDER,
        HAS_ANALYSIS
    }
    public void startDatabase(){
        log.log(Level.INFO, "Starting database ...");

        graphDb = new GraphDatabaseFactory()
                .newEmbeddedDatabaseBuilder(neo4jDBPath)
                .newGraphDatabase();

        Neo4j.registerShutdownHook(graphDb);
    }
    public void createIndexes() {
        log.log(Level.INFO, "Adding constraints ...");

        Neo4j.createConstraint(graphDb, sampleLabel, "SampleId");
        Neo4j.createConstraint(graphDb, variantLabel, "VariantId");
        Neo4j.createConstraint(graphDb, variantLabel, "RsId");
        Neo4j.createConstraint(graphDb, featureLabel, "FeatureId");
        Neo4j.createConstraint(graphDb, canonicalLabel, "FeatureId");
        Neo4j.createConstraint(graphDb, symbolLabel, "SymbolId");
        Neo4j.createConstraint(graphDb, disorderLabel, "Title");
        Neo4j.createConstraint(graphDb, runInfoLabel, "AnalysisId");
        Neo4j.createIndex(graphDb, runInfoLabel, "LibraryId");
        Neo4j.createIndex(graphDb, runInfoLabel, "RunId");
    }

    public void loadVCFFiles() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Loading VCFs into memory ...");

        int n = 0;
        Iterator<VariantContext> variantContextIterator;

        //read variant VCF file
        variantContextIterator = variantVcfFileReader.iterator();
        while (variantContextIterator.hasNext()) {
            VariantContext variant = variantContextIterator.next();
            ++n;

            if (n > 250){
                break;
            }

            if (!variant.isFiltered() && variant.isVariant()){
                vcfBody.add(variant);
            }

        }

        //read annotation VCF file
        variantContextIterator = annotationVcfFileReader.iterator();
        while (variantContextIterator.hasNext()){
            VariantContext variant = variantContextIterator.next();

            if (variant.getAlternateAlleles().size() != 1) throw new InvalidPropertiesFormatException("Allele " + variant.getAlleles().toString() + " is not diploid");

            GenomeVariant genomeVariant = new GenomeVariant(variant.getContig(), variant.getStart(), variant.getReference().getBaseString(), variant.getAlternateAllele(0).getBaseString());

            if (!vepAnnotations.containsKey(genomeVariant)) vepAnnotations.put(genomeVariant, new HashSet<VEPAnnotation>());

            //split annotations and make unique
            try {

                //one annotation
                VEPAnnotation vepAnnotation = new VEPAnnotation((String) variant.getAttribute("CSQ"));
                vepAnnotation.parseAnnotation();

                vepAnnotations.get(genomeVariant).add(vepAnnotation);

            } catch (ClassCastException e) {

                //multiple annotations
                for (String annotation : (ArrayList<String>) variant.getAttribute("CSQ")) {

                    VEPAnnotation vepAnnotation = new VEPAnnotation(annotation);
                    vepAnnotation.parseAnnotation();

                    vepAnnotations.get(genomeVariant).add(vepAnnotation);

                }
            }

        }

    }
    public void addSampleAndRunInfoNodes() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Adding sample and run info nodes ...");

        HashMap<String, Object> properties = new HashMap<>();
        ArrayList<String> sampleIds = variantVcfFileReader.getFileHeader().getSampleNamesInOrder();

        String libraryId = "K15-0000"; //TODO extract from VCF - Shire worklist
        String runId = "150716_D00501_0047_BHB092ADXX"; //TODO extract from VCF - flowcell/chipId
        //TODO get worksheet position (n)

        for (int n = 0; n < sampleIds.size(); ++n){

            //add sample
            Node sampleIdNode = Neo4j.matchOrCreateUniqueNode(graphDb, sampleLabel, "SampleId", sampleIds.get(n));

            //add run info
            properties.put("LibraryId", libraryId);
            properties.put("RunId", runId);
            properties.put("SampleNo", n);
            properties.put("AnalysisId", libraryId + "_" + n + "_" + runId);

            Node runInfoNode = Neo4j.addNode(graphDb, runInfoLabel, properties);
            properties.clear();

            //link sample and runInfo
            Neo4j.createRelationship(graphDb, sampleIdNode, runInfoNode, relTypes.HAS_ANALYSIS, properties, false);
            runInfoNodes.put(sampleIds.get(n), runInfoNode);
        }

    }
    public void addVariantNodesAndGenotypeRelationships() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Adding variants and genotypes ...");

        GenomeVariant genomeVariant;

        //read VCF records
        for (VariantContext vcfRecord : vcfBody) {
            Iterator<Genotype> genotypeIterator = vcfRecord.getGenotypes().iterator();

            while (genotypeIterator.hasNext()) {
                Genotype genotype = genotypeIterator.next();

                if (genotype.isNoCall() || genotype.isHomRef() || genotype.isMixed()) continue;
                if (genotype.getPloidy() != 2 || genotype.getAlleles().size() != 2) throw new InvalidPropertiesFormatException("Allele " + genotype.getAlleles().toString() + " is not diploid");

                //add new variants to the DB
                if (genotype.isHom()){

                    if (genotype.getAlleles().get(1).getBaseString().equals("*")) continue;

                    genomeVariant = new GenomeVariant(vcfRecord.getContig(), vcfRecord.getStart(), vcfRecord.getReference().getBaseString(), genotype.getAlleles().get(1).getBaseString());
                    genomeVariant.convertToMinimalRepresentation();

                    addVariantNodesAndGenotypeRelationshipsHelper(genomeVariant, genotype.getGQ(), runInfoNodes.get(genotype.getSampleName()), relTypes.HAS_HOM_VARIANT);

                } else if (genotype.isHet()){

                    if (genotype.getAlleles().get(1).getBaseString().equals("*")) continue;

                    genomeVariant = new GenomeVariant(vcfRecord.getContig(), vcfRecord.getStart(), vcfRecord.getReference().getBaseString(), genotype.getAlleles().get(1).getBaseString());
                    genomeVariant.convertToMinimalRepresentation();

                    addVariantNodesAndGenotypeRelationshipsHelper(genomeVariant, genotype.getGQ(), runInfoNodes.get(genotype.getSampleName()), relTypes.HAS_HET_VARIANT);

                    if (genotype.isHetNonRef()){

                        if (genotype.getAlleles().get(0).getBaseString().equals("*")) continue;

                        genomeVariant = new GenomeVariant(vcfRecord.getContig(), vcfRecord.getStart(), vcfRecord.getReference().getBaseString(), genotype.getAlleles().get(0).getBaseString());
                        genomeVariant.convertToMinimalRepresentation();

                        addVariantNodesAndGenotypeRelationshipsHelper(genomeVariant, genotype.getGQ(), runInfoNodes.get(genotype.getSampleName()), relTypes.HAS_HET_VARIANT);

                    }

                } else {
                    throw new InvalidPropertiesFormatException("Dosage unknown: " + vcfRecord.toString());
                }

            }

        }

    }
    public void addPopulationFrequencies(){
        log.log(Level.INFO, "Adding population frequencies ...");

        //TODO and dbSnp RsID
    }
    public void addAnnotations() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Adding annotations ...");

        HashMap<String, Object> properties = new HashMap<>();

        for (Map.Entry<GenomeVariant, Node> variant : newVariantNodes.entrySet()) { //loop over new variants
            for (VEPAnnotation annotation : vepAnnotations.get(variant.getKey())){

                if (!annotationNodes.containsKey(variant.getKey())) annotationNodes.put(variant.getKey(), new HashMap<String, Node>());

                //add symbols
                if (annotation.getSymbol() != null && !annotation.getSymbol().equals("") && !symbolNodes.containsKey(annotation.getSymbol())){
                    ArrayList<Node> nodes = Neo4j.getNodes(graphDb, symbolLabel, "SymbolId", annotation.getSymbol());

                    if (nodes.size() == 0){

                        //add symbol
                        properties.put("SymbolId", annotation.getSymbol());
                        if (annotation.getSymbolSource() != null && !annotation.getSymbolSource().equals("")) properties.put("SymbolSource", annotation.getSymbolSource());

                        symbolNodes.put(annotation.getSymbol(), Neo4j.addNode(graphDb, symbolLabel, properties));

                    } else {
                        symbolNodes.put(annotation.getSymbol(), nodes.get(0));
                    }

                }

                //add feature
                if (annotation.getFeature() != null && !annotation.getFeature().equals("") && !featureNodes.containsKey(annotation.getFeature())){
                    ArrayList<Node> nodes = Neo4j.getNodes(graphDb, featureLabel, "FeatureId", annotation.getFeature());

                    if (nodes.size() == 0) {

                        properties.clear();
                        if(annotation.getFeatureType() != null) properties.put("FeatureType", annotation.getFeatureType());
                        if(annotation.getBiotype() != null) properties.put("Biotype", annotation.getBiotype());
                        if(annotation.getStrand() != null) properties.put("Strand", annotation.getStrand());
                        if(annotation.getExon() != null) properties.put("TotalExons", annotation.getExon().split("/")[1]);

                        featureNodes.put(annotation.getFeature(), Neo4j.addNode(graphDb, featureLabel, properties));

                        if (annotation.getCanonical() != null && annotation.getCanonical().equals("YES")) Neo4j.addNodeLabel(graphDb, featureNodes.get(annotation.getFeature()), canonicalLabel);

                    } else {
                        featureNodes.put(annotation.getFeature(), nodes.get(0));
                    }

                }

                //add functional annotations
                properties.clear();
                if(annotation.getHgvsCoding() != null) properties.put("HGVSc", annotation.getHgvsCoding());
                if(annotation.getHgvsProtein() != null) properties.put("HGVSp", annotation.getHgvsProtein());
                if(annotation.getExon() != null) properties.put("Exon", annotation.getExon().split("/")[0]);
                if(annotation.getIntron() != null) properties.put("Intron", annotation.getIntron().split("/")[0]);

                annotationNodes.get(variant.getKey()).put(annotation.getFeature(), Neo4j.addNode(graphDb, annotationLabel, properties));

                //add omim disorders and link to disease
                properties.clear();
                if (geneMap2.containsKey(annotation.getSymbol())) {
                    for (String disorderTitle : geneMap2.get(annotation.getSymbol())){
                        Node node = Neo4j.matchOrCreateUniqueNode(graphDb, disorderLabel, "Title", disorderTitle);
                        Neo4j.createRelationship(graphDb, symbolNodes.get(annotation.getSymbol()), node, relTypes.HAS_ASSOCIATED_DISORDER, properties, false);
                    }
                }

                //link consequences
                properties.clear();
                if (annotation.getConsequences().size() > 0){
                    for (String consequence : annotation.getConsequences()){
                        Neo4j.createRelationship(graphDb, newVariantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), DynamicRelationshipType.withName("HAS_" + consequence.toUpperCase() + "_CONSEQUENCE"), properties, false);
                    }
                } else {
                    Neo4j.createRelationship(graphDb, newVariantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_UNKNOWN_CONSEQUENCE, properties, false);
                }

                //add in feature relationship
                if (annotationNodes.get(variant.getKey()).get(annotation.getFeature()) != null && featureNodes.get(annotation.getFeature()) != null){
                    Neo4j.createRelationship(graphDb, annotationNodes.get(variant.getKey()).get(annotation.getFeature()), featureNodes.get(annotation.getFeature()), relTypes.IN_FEATURE, new HashMap<String, Object>(), false);
                }

                //add in symbol relationship
                if (featureNodes.get(annotation.getFeature()) != null && symbolNodes.get(annotation.getSymbol()) != null){
                    Neo4j.createRelationship(graphDb, featureNodes.get(annotation.getFeature()), symbolNodes.get(annotation.getSymbol()), relTypes.IN_SYMBOL, new HashMap<String, Object>(), false);
                }

            }
        }

    }

    private void addVariantNodesAndGenotypeRelationshipsHelper(GenomeVariant genomeVariant, int genotypeQuality, Node runInfoNode, RelationshipType relationshipType){

        HashMap<String, Object> properties = new HashMap<>();

        if (!newVariantNodes.containsKey(genomeVariant) && !existingNodes.containsKey(genomeVariant)) {
            ArrayList<Node> nodes = Neo4j.getNodes(graphDb, variantLabel, "VariantId", genomeVariant.getConcatenatedVariant());

            if (nodes.size() == 0){
                properties.put("VariantId", genomeVariant.getConcatenatedVariant());
                newVariantNodes.put(genomeVariant, Neo4j.addNode(graphDb, variantLabel, properties));
            } else {
                existingNodes.put(genomeVariant, nodes.get(0));
            }
        }
        properties.clear();

        //create genotype relationship
        properties.put("Quality", genotypeQuality);
        Neo4j.createRelationship(graphDb, runInfoNode, newVariantNodes.get(genomeVariant), relationshipType, properties, false);

    }

    public void shutdownDatabase(){
        log.log(Level.INFO, "Shutting down database ...");
        graphDb.shutdown();
    }

}