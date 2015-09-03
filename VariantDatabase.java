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

//NB - assumption = if variant exists in DB so do the annotations
//TODO convert VCF to minimal represented and remove sample-specific info; annotate then use as input
//TODO check reference genome for 1kgp3 against b37
//TODO check all variants have annotations
//TODO import dbsnfp and use e v79
//TODO import phenotype associated with gene
//TODO retrieve pubmed abstracts (web ui)
//TODO check alamut for extra functionality
//TODO skip NTC
//TODO add variant type label

public class VariantDatabase {

    private static final Logger log = Logger.getLogger(VariantDatabase.class.getName());

    private File neo4jDBPath;
    private GraphDatabaseService graphDb;
    private VCFFileReader variantVcfFileReader, annotationVcfFileReader;
    private ArrayList<VariantContext> vcfBody = new ArrayList<>(); //VCF file body
    private HashMap<GenomeVariant, HashSet<VEPAnnotation>> vepAnnotations = new HashMap<>(); //all VEP annotations
    private HashMap<GenomeVariant, HashMap<String, Double>> populationFrequencies = new HashMap<>(); //population frequencies from mulitple sources
    private HashMap<GenomeVariant, Node> newVariantNodes = new HashMap<>(); //new variants added during this session
    private HashMap<GenomeVariant, Node> existingNodes = new HashMap<>();
    private HashMap<String, Node> sampleNodes = new HashMap<>(); //samples added during this session
    private HashMap<GenomeVariant, HashMap<String, Node>> annotationNodes = new HashMap<>(); //
    private HashMap<String, Node> symbolNodes = new HashMap<>(); //symbols added during this session
    private HashMap<String, Node> featureNodes = new HashMap<>(); //features added during this session

    //TODO add variables to VCF
    private String libraryId = "K15-0000"; //TODO extract from VCF - Shire worklist
    private String runId = "150716_D00501_0047_BHB092ADXX"; //TODO extract from VCF - flowcell/chipId
    private int sampleNo = 0; //TODO extract from VCF - index number/position on sampleSheet

    private Label sampleLabel = DynamicLabel.label("Sample");
    private Label variantLabel = DynamicLabel.label("Variant");
    private Label annotationLabel = DynamicLabel.label("Annotation");
    private Label symbolLabel = DynamicLabel.label("Symbol");
    private Label canonicalLabel = DynamicLabel.label("Canonical");
    private Label featureLabel = DynamicLabel.label("Feature");

    public VariantDatabase(VCFFileReader variantVcfFileReader, VCFFileReader annotationVcfFileReader, File neo4jDBPath){
        this.variantVcfFileReader = variantVcfFileReader;
        this.neo4jDBPath = neo4jDBPath;
        this.annotationVcfFileReader = annotationVcfFileReader;
    }

    private enum relTypes implements RelationshipType
    {
        HAS_HET_VARIANT,
        HAS_HOM_VARIANT,
        IN_SYMBOL,
        IN_FEATURE
    }
    public void startDatabase(){
        log.log(Level.INFO, "Starting database ...");

        //create DB
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

            if (n > 100){
                //break;
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
    public void addSampleNodes() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Adding sample nodes ...");

        //TODO add sample type

        for (String sampleId : variantVcfFileReader.getFileHeader().getSampleNamesInOrder()){
            sampleNodes.put(sampleId, Neo4j.matchOrCreateUniqueNode(graphDb, sampleLabel, "SampleId", sampleId));
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

                    addVariantNodesAndGenotypeRelationshipsHelper(genomeVariant, genotype.getGQ(), sampleNodes.get(genotype.getSampleName()), relTypes.HAS_HOM_VARIANT);

                } else if (genotype.isHet()){

                    if (genotype.getAlleles().get(1).getBaseString().equals("*")) continue;

                    genomeVariant = new GenomeVariant(vcfRecord.getContig(), vcfRecord.getStart(), vcfRecord.getReference().getBaseString(), genotype.getAlleles().get(1).getBaseString());
                    genomeVariant.convertToMinimalRepresentation();

                    addVariantNodesAndGenotypeRelationshipsHelper(genomeVariant, genotype.getGQ(), sampleNodes.get(genotype.getSampleName()), relTypes.HAS_HET_VARIANT);

                    if (genotype.isHetNonRef()){

                        if (genotype.getAlleles().get(0).getBaseString().equals("*")) continue;

                        genomeVariant = new GenomeVariant(vcfRecord.getContig(), vcfRecord.getStart(), vcfRecord.getReference().getBaseString(), genotype.getAlleles().get(0).getBaseString());
                        genomeVariant.convertToMinimalRepresentation();

                        addVariantNodesAndGenotypeRelationshipsHelper(genomeVariant, genotype.getGQ(), sampleNodes.get(genotype.getSampleName()), relTypes.HAS_HET_VARIANT);

                    }

                } else {
                    throw new InvalidPropertiesFormatException("Dosage unknown: " + vcfRecord.toString());
                }

            }

        }

    }
    public void addPopulationFrequencies(){
        //TODO and dbSnp RsID
    }
    public void addSymbolNodes() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Adding symbols ...");

        //loop over new variants
        for (Map.Entry<GenomeVariant, Node> variant : newVariantNodes.entrySet()) {
            if (vepAnnotations.containsKey(variant.getKey())) {

                for (VEPAnnotation annotation : vepAnnotations.get(variant.getKey())){

                    //skip symbols already imported during this session
                    if (annotation.getSymbol() != null && !annotation.getSymbol().equals("") && !symbolNodes.containsKey(annotation.getSymbol())){
                        symbolNodes.put(annotation.getSymbol(), Neo4j.matchOrCreateUniqueNode(graphDb, symbolLabel, "SymbolId", annotation.getSymbol()));
                    }

                }

            }
        }

    }
    public void addFeatureNodes() {
        log.log(Level.INFO, "Adding features ...");

        //loop over new variants
        for (Map.Entry<GenomeVariant, Node> variant : newVariantNodes.entrySet()) {
            if (vepAnnotations.containsKey(variant.getKey())) {

                for (VEPAnnotation annotation : vepAnnotations.get(variant.getKey())){

                    if (annotation.getFeature() != null && !annotation.getFeature().equals("") && !featureNodes.containsKey(annotation.getFeature())){ //skip features already imported during this session

                        ArrayList<Node> nodes = Neo4j.getNodes(graphDb, featureLabel, "FeatureId", annotation.getFeature());

                        if (nodes.size() == 0) {

                            HashMap<String, Object> properties = new HashMap<>();

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

                }

            }
        }

    }
    public void addFunctionalAnnotationNodes() {
        log.log(Level.INFO, "Adding functional annotations ...");

        //TODO add more annotations from dbsnfp

        //loop over variants
        for (Map.Entry<GenomeVariant, Node> variant : newVariantNodes.entrySet()){
            if (vepAnnotations.containsKey(variant.getKey())){

                if (!annotationNodes.containsKey(variant.getKey())) annotationNodes.put(variant.getKey(), new HashMap<String, Node>());

                //loop over functional annotations for this variant
                for (VEPAnnotation annotation : vepAnnotations.get(variant.getKey())) {

                    HashMap<String, Object> properties = new HashMap<>();

                    if(annotation.getHgvsCoding() != null) properties.put("HGVSc", annotation.getHgvsCoding());
                    if(annotation.getHgvsProtein() != null) properties.put("HGVSp", annotation.getHgvsProtein());
                    if(annotation.getPolyphen() != null) properties.put("PolyPhen", annotation.getPolyphen());
                    if(annotation.getSift() != null) properties.put("SIFT", annotation.getSift());
                    if(annotation.getExon() != null) properties.put("Exon", annotation.getExon().split("/")[0]);
                    if(annotation.getIntron() != null) properties.put("Intron", annotation.getIntron().split("/")[0]);

                    annotationNodes.get(variant.getKey()).put(annotation.getFeature(), Neo4j.addNode(graphDb, annotationLabel, properties));
                }

            }
        }

    }
    public void addConsequenceRelationships() {
        log.log(Level.INFO, "Linking variants To functional annotations ...");

        //loop over variants
        for (Map.Entry<GenomeVariant, Node> variant : newVariantNodes.entrySet()){
            if (vepAnnotations.containsKey(variant.getKey())){

                HashMap<String, Object> properties = new HashMap<>();

                //loop over all annotations by for this allele
                for (VEPAnnotation annotation : vepAnnotations.get(variant.getKey())){

                    for (String consequence : annotation.getConsequences()){

                        //create relationship type
                        Neo4j.createRelationship(graphDb, newVariantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), DynamicRelationshipType.withName("HAS_" + consequence.toUpperCase() + "_CONSEQUENCE"), properties, false);
                    }

                }
            }

        }
    }
    public void addInFeatureRelationships() {
        log.log(Level.INFO, "Linking annotations to features ...");

        //loop over variants
        for (Map.Entry<GenomeVariant, Node> variant : newVariantNodes.entrySet()){
            if (vepAnnotations.containsKey(variant.getKey())){

                //loop over functional annotations for this variant
                for (VEPAnnotation annotation : vepAnnotations.get(variant.getKey())) {
                    Neo4j.createRelationship(graphDb, annotationNodes.get(variant.getKey()).get(annotation.getFeature()), featureNodes.get(annotation.getFeature()), relTypes.IN_FEATURE, new HashMap<String, Object>(), false);
                }

            }
        }

    }
    public void addInSymbolRelationships() {
        log.log(Level.INFO, "Linking features to symbol ...");

        //loop over variants
        for (Map.Entry<GenomeVariant, Node> variant : newVariantNodes.entrySet()) {
            if (vepAnnotations.containsKey(variant.getKey())) {

                //loop over functional annotations for this variant
                for (VEPAnnotation annotation : vepAnnotations.get(variant.getKey())) {
                    Neo4j.createRelationship(graphDb, featureNodes.get(annotation.getFeature()), symbolNodes.get(annotation.getSymbol()), relTypes.IN_SYMBOL, new HashMap<String, Object>(), false);
                }

            }
        }

    }

    private void addVariantNodesAndGenotypeRelationshipsHelper(GenomeVariant genomeVariant, int genotypeQuality, Node sampleNode, RelationshipType relationshipType){

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
        properties.put("LibraryId", libraryId);
        properties.put("RunId", runId);
        properties.put("SampleNo", sampleNo);
        properties.put("Quality", genotypeQuality);

        if (!Neo4j.isNeighbourNodeWithSuppliedProperties(graphDb, sampleNode, newVariantNodes.get(genomeVariant), Direction.OUTGOING, relationshipType, properties)){
            Neo4j.createRelationship(graphDb, sampleNode, newVariantNodes.get(genomeVariant), relationshipType, properties, false);
        }

    }

    public void shutdownDatabase(){
        log.log(Level.INFO, "Shutting down database ...");
        graphDb.shutdown();
    }

}