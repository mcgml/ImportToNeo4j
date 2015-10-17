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

//TODO retrieve pubmed abstracts (web ui)
//TODO check alamut for extra functionality
//TODO add sample type
//TODO add pipeline name & version
//todo add ddd afs
//todo add clinvar
//todo deal with mixed genotypes
//todo add shire info and seq info

public class VariantDatabase {
    private static final Logger log = Logger.getLogger(VariantDatabase.class.getName());

    private File neo4jDBPath;
    private GraphDatabaseService graphDb;
    private VCFFileReader variantVcfFileReader, annotationVcfFileReader;
    private ArrayList<VariantContext> vcfBody = new ArrayList<>(); //VCF file body
    private HashMap<GenomeVariant, HashSet<VEPAnnotationv75>> vepAnnotations = new HashMap<>(); //all VEP annotations
    private HashMap<GenomeVariant, HashMap<String, Object>> genomeVariantAnnotations = new HashMap<>(); //population frequencies from multiple sources
    private HashMap<GenomeVariant, Node> newVariantNodes = new HashMap<>(); //new variants added during this session
    private HashMap<GenomeVariant, Node> existingNodes = new HashMap<>();
    private HashMap<String, Node> runInfoNodes = new HashMap<>(); //samples added during this session
    private HashMap<GenomeVariant, HashMap<String, Node>> annotationNodes = new HashMap<>(); //variant = Transcript:annotationNode
    private HashMap<String, Node> symbolNodes = new HashMap<>(); //symbols added during this session
    private HashMap<String, Node> featureNodes = new HashMap<>(); //features added during this session
    private Node userNode;

    private String libraryId = "K15-0000"; //TODO extract from VCF - Shire worklist
    private String runId = "150716_D00501_0047_BHB092ADXX"; //TODO extract from VCF - flowcell/chipId
    private String panelName = "Illumina TruSight One"; //TODO extract from VCF - assay type
    private String pipelineName = "IlluminaTruSightOne"; //todo get from VCF
    private int pipelineVersion = 1; //todo get from VCF

    private int minimumAlellesForAFCalculation = 120;

    public VariantDatabase(VCFFileReader variantVcfFileReader, VCFFileReader annotationVcfFileReader, File neo4jDBPath){
        this.variantVcfFileReader = variantVcfFileReader;
        this.neo4jDBPath = neo4jDBPath;
        this.annotationVcfFileReader = annotationVcfFileReader;
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

        Neo4j.createConstraint(graphDb, Neo4j.getSampleLabel(), "SampleId");
        Neo4j.createIndex(graphDb, Neo4j.getRunInfoLabel(), "LibraryId");
        Neo4j.createIndex(graphDb, Neo4j.getRunInfoLabel(), "RunId");
        Neo4j.createConstraint(graphDb, Neo4j.getRunInfoLabel(), "AnalysisId");
        Neo4j.createConstraint(graphDb, Neo4j.getVariantLabel(), "VariantId");
        Neo4j.createConstraint(graphDb, Neo4j.getMtChromosomeLabel(), "VariantId");
        Neo4j.createConstraint(graphDb, Neo4j.getXChromosomeLabel(), "VariantId");
        Neo4j.createConstraint(graphDb, Neo4j.getYChromosomeLabel(), "VariantId");
        Neo4j.createConstraint(graphDb, Neo4j.getAutoChromosomeLabel(), "VariantId");
        Neo4j.createIndex(graphDb, Neo4j.getVariantLabel(), "Id");
        Neo4j.createConstraint(graphDb, Neo4j.getFeatureLabel(), "FeatureId");
        Neo4j.createConstraint(graphDb, Neo4j.getCanonicalLabel(), "FeatureId");
        Neo4j.createConstraint(graphDb, Neo4j.getSymbolLabel(), "SymbolId");
        Neo4j.createConstraint(graphDb, Neo4j.getVirtualPanelLabel(), "PanelName");
        Neo4j.createConstraint(graphDb, Neo4j.getUserLabel(), "UserName");

    }
    public void addUsers() throws InvalidPropertiesFormatException {
        userNode = Neo4j.matchOrCreateUniqueNode(graphDb, Neo4j.getUserLabel(), "UserName", "ml");
    }

    public void loadVCFFiles() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Loading VCFs into memory ...");

        //read variant VCF file
        Iterator<VariantContext> variantContextIterator = variantVcfFileReader.iterator();

        while (variantContextIterator.hasNext()) {
            VariantContext variant = variantContextIterator.next();

            if (!variant.isFiltered() && variant.isVariant()){
                vcfBody.add(variant);
            }

        }

        //read annotation VCF file
        Iterator<VariantContext> annotationContextIterator = annotationVcfFileReader.iterator();

        while (annotationContextIterator.hasNext()){
            VariantContext variant = annotationContextIterator.next();

            if (variant.getAlternateAlleles().size() != 1) throw new InvalidPropertiesFormatException("Allele " + variant.getAlleles().toString() + " is not diploid");

            GenomeVariant genomeVariant = new GenomeVariant(variant.getContig(), variant.getStart(), variant.getReference().getBaseString(), variant.getAlternateAllele(0).getBaseString());

            if (!vepAnnotations.containsKey(genomeVariant)) vepAnnotations.put(genomeVariant, new HashSet<VEPAnnotationv75>());
            if (!genomeVariantAnnotations.containsKey(genomeVariant)) genomeVariantAnnotations.put(genomeVariant, new HashMap<String, Object>());

            //split annotations and make unique
            try {

                //one annotation
                VEPAnnotationv75 vepAnnotationv75 = new VEPAnnotationv75((String) variant.getAttribute("CSQ"));
                vepAnnotationv75.parseAnnotation();

                vepAnnotations.get(genomeVariant).add(vepAnnotationv75);

            } catch (ClassCastException e) {

                //multiple annotations
                for (String annotation : (ArrayList<String>) variant.getAttribute("CSQ")) {

                    VEPAnnotationv75 vepAnnotationv75 = new VEPAnnotationv75(annotation);
                    vepAnnotationv75.parseAnnotation();

                    vepAnnotations.get(genomeVariant).add(vepAnnotationv75);

                }
            }

            //add population frequencies
            if (variant.getAttribute("1kgp3.EAS_AF") != null && !variant.getAttribute("1kgp3.EAS_AF").equals(".")) {
                genomeVariantAnnotations.get(genomeVariant).put("onekGPhase3_EAS_AF", Float.parseFloat((String) variant.getAttribute("1kgp3.EAS_AF")));
            }
            if (variant.getAttribute("1kgp3.EUR_AF") != null && !variant.getAttribute("1kgp3.EUR_AF").equals(".")) {
                genomeVariantAnnotations.get(genomeVariant).put("onekGPhase3_EUR_AF", Float.parseFloat((String) variant.getAttribute("1kgp3.EUR_AF")));
            }
            if (variant.getAttribute("1kgp3.AFR_AF") != null && !variant.getAttribute("1kgp3.AFR_AF").equals(".")) {
                genomeVariantAnnotations.get(genomeVariant).put("onekGPhase3_AFR_AF", Float.parseFloat((String) variant.getAttribute("1kgp3.AFR_AF")));
            }
            if (variant.getAttribute("1kgp3.AMR_AF") != null && !variant.getAttribute("1kgp3.AMR_AF").equals(".")) {
                genomeVariantAnnotations.get(genomeVariant).put("onekGPhase3_AMR_AF", Float.parseFloat((String) variant.getAttribute("1kgp3.AMR_AF")));
            }
            if (variant.getAttribute("1kgp3.SAS_AF") != null && !variant.getAttribute("1kgp3.SAS_AF").equals(".")) {
                genomeVariantAnnotations.get(genomeVariant).put("onekGPhase3_SAS_AF", Float.parseFloat((String) variant.getAttribute("1kgp3.SAS_AF")));
            }
            if (variant.getAttribute("ExAC.AC_AFR") != null && !variant.getAttribute("ExAC.AC_AFR").equals(".") && variant.getAttribute("ExAC.AN_AFR") != null && !variant.getAttribute("ExAC.AN_AFR").equals(".") && Integer.parseInt((String) variant.getAttribute("ExAC.AN_AFR")) > minimumAlellesForAFCalculation) {
                genomeVariantAnnotations.get(genomeVariant).put("ExAC_AFR_AF", Float.parseFloat((String) variant.getAttribute("ExAC.AC_AFR")) / Float.parseFloat((String) variant.getAttribute("ExAC.AN_AFR")));
            }
            if (variant.getAttribute("ExAC.AC_AMR") != null && !variant.getAttribute("ExAC.AC_AMR").equals(".") && variant.getAttribute("ExAC.AN_AMR") != null && !variant.getAttribute("ExAC.AN_AMR").equals(".") && Integer.parseInt((String) variant.getAttribute("ExAC.AN_AMR")) > minimumAlellesForAFCalculation) {
                genomeVariantAnnotations.get(genomeVariant).put("ExAC_AMR_AF", Float.parseFloat((String) variant.getAttribute("ExAC.AC_AMR")) / Float.parseFloat((String) variant.getAttribute("ExAC.AN_AMR")));
            }
            if (variant.getAttribute("ExAC.AC_EAS") != null && !variant.getAttribute("ExAC.AC_EAS").equals(".") && variant.getAttribute("ExAC.AN_EAS") != null && !variant.getAttribute("ExAC.AN_EAS").equals(".") && Integer.parseInt((String) variant.getAttribute("ExAC.AN_EAS")) > minimumAlellesForAFCalculation) {
                genomeVariantAnnotations.get(genomeVariant).put("ExAC_EAS_AF", Float.parseFloat((String) variant.getAttribute("ExAC.AC_EAS")) / Float.parseFloat((String) variant.getAttribute("ExAC.AN_EAS")));
            }
            if (variant.getAttribute("ExAC.AC_FIN") != null && !variant.getAttribute("ExAC.AC_FIN").equals(".") && variant.getAttribute("ExAC.AN_FIN") != null && !variant.getAttribute("ExAC.AN_FIN").equals(".") && Integer.parseInt((String) variant.getAttribute("ExAC.AN_FIN")) > minimumAlellesForAFCalculation) {
                genomeVariantAnnotations.get(genomeVariant).put("ExAC_FIN_AF", Float.parseFloat((String) variant.getAttribute("ExAC.AC_FIN")) / Float.parseFloat((String) variant.getAttribute("ExAC.AN_FIN")));
            }
            if (variant.getAttribute("ExAC.AC_NFE") != null && !variant.getAttribute("ExAC.AC_NFE").equals(".") && variant.getAttribute("ExAC.AN_NFE") != null && !variant.getAttribute("ExAC.AN_NFE").equals(".") && Integer.parseInt((String) variant.getAttribute("ExAC.AN_NFE")) > minimumAlellesForAFCalculation) {
                genomeVariantAnnotations.get(genomeVariant).put("ExAC_NFE_AF", Float.parseFloat((String) variant.getAttribute("ExAC.AC_NFE")) / Float.parseFloat((String) variant.getAttribute("ExAC.AN_NFE")));
            }
            if (variant.getAttribute("ExAC.AC_OTH") != null && !variant.getAttribute("ExAC.AC_OTH").equals(".") && variant.getAttribute("ExAC.AN_OTH") != null && !variant.getAttribute("ExAC.AN_OTH").equals(".") && Integer.parseInt((String) variant.getAttribute("ExAC.AN_OTH")) > minimumAlellesForAFCalculation) {
                genomeVariantAnnotations.get(genomeVariant).put("ExAC_OTH_AF", Float.parseFloat((String) variant.getAttribute("ExAC.AC_OTH")) / Float.parseFloat((String) variant.getAttribute("ExAC.AN_OTH")));
            }
            if (variant.getAttribute("ExAC.AC_SAS") != null && !variant.getAttribute("ExAC.AC_SAS").equals(".") && variant.getAttribute("ExAC.AN_SAS") != null && !variant.getAttribute("ExAC.AN_SAS").equals(".") && Integer.parseInt((String) variant.getAttribute("ExAC.AN_SAS")) > minimumAlellesForAFCalculation) {
                genomeVariantAnnotations.get(genomeVariant).put("ExAC_SAS_AF", Float.parseFloat((String) variant.getAttribute("ExAC.AC_SAS")) / Float.parseFloat((String) variant.getAttribute("ExAC.AN_SAS")));
            }

            //add dbSNP id
            if (!variant.getID().equals(".") && !variant.getID().equals("")){
                genomeVariantAnnotations.get(genomeVariant).put("Id", variant.getID());
            }

        }

    }
    public void addSampleAndRunInfoNodes() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Adding sample and run info nodes ...");

        HashMap<String, Object> properties = new HashMap<>();
        ArrayList<String> sampleIds = variantVcfFileReader.getFileHeader().getSampleNamesInOrder();

        for (short n = 0; n < sampleIds.size(); ++n){

            //add sample
            Node sampleIdNode = Neo4j.matchOrCreateUniqueNode(graphDb, Neo4j.getSampleLabel(), "SampleId", sampleIds.get(n));

            //add run info
            properties.put("LibraryId", libraryId);
            properties.put("RunId", runId);
            properties.put("SampleNo", n);
            properties.put("AnalysisId", libraryId + "_" + n + "_" + runId);
            properties.put("PanelName", panelName);
            properties.put("PipelineName", pipelineName);
            properties.put("PipelineVersion", pipelineVersion);

            Node runInfoNode = Neo4j.addNode(graphDb, Neo4j.getRunInfoLabel(), properties); //need to crash if already exists
            properties.clear();

            //link sample and runInfo
            Neo4j.createRelationship(graphDb, sampleIdNode, runInfoNode, Neo4j.getHasAnalysisRelationship(), properties, false);
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

                if (genotype.isNoCall() || genotype.isHomRef()) continue;
                if (genotype.isMixed()){
                    log.log(Level.WARNING, genotype.getSampleName() + ": " + vcfRecord.getContig() + " " + vcfRecord.getStart() + " " + vcfRecord.getReference() + vcfRecord.getAlternateAlleles().toString() + " has mixed genotype and could not be added.");
                    continue;
                }
                if (genotype.getPloidy() != 2 || genotype.getAlleles().size() != 2) throw new InvalidPropertiesFormatException("Allele " + genotype.getAlleles().toString() + " is not diploid");

                //add new variants to the DB
                if (genotype.isHom()){

                    if (genotype.getAlleles().get(1).getBaseString().equals("*")) continue;

                    genomeVariant = new GenomeVariant(vcfRecord.getContig(), vcfRecord.getStart(), vcfRecord.getReference().getBaseString(), genotype.getAlleles().get(1).getBaseString());
                    genomeVariant.convertToMinimalRepresentation();

                    addVariantNodesAndGenotypeRelationshipsHelper(genomeVariant, (short) genotype.getGQ(), runInfoNodes.get(genotype.getSampleName()), Neo4j.getHasHomVariantRelationship());

                } else if (genotype.isHet()){

                    if (genotype.getAlleles().get(1).getBaseString().equals("*")) continue;

                    genomeVariant = new GenomeVariant(vcfRecord.getContig(), vcfRecord.getStart(), vcfRecord.getReference().getBaseString(), genotype.getAlleles().get(1).getBaseString());
                    genomeVariant.convertToMinimalRepresentation();

                    addVariantNodesAndGenotypeRelationshipsHelper(genomeVariant, (short) genotype.getGQ(), runInfoNodes.get(genotype.getSampleName()), Neo4j.getHasHetVariantRelationship());

                    if (genotype.isHetNonRef()){

                        if (genotype.getAlleles().get(0).getBaseString().equals("*")) continue;

                        genomeVariant = new GenomeVariant(vcfRecord.getContig(), vcfRecord.getStart(), vcfRecord.getReference().getBaseString(), genotype.getAlleles().get(0).getBaseString());
                        genomeVariant.convertToMinimalRepresentation();

                        addVariantNodesAndGenotypeRelationshipsHelper(genomeVariant, (short) genotype.getGQ(), runInfoNodes.get(genotype.getSampleName()), Neo4j.getHasHetVariantRelationship());

                    }

                } else {
                    throw new InvalidPropertiesFormatException("Inheritance unknown: " + vcfRecord.toString());
                }

            }

        }

    }
    public void addAnnotations() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Adding annotations ...");

        HashMap<String, Object> properties = new HashMap<>();

        //loop over new variants adding during this session and add annotations
        for (Map.Entry<GenomeVariant, Node> variant : newVariantNodes.entrySet()) {

            //check variant has annotation
            if (!vepAnnotations.containsKey(variant.getKey())){
                log.log(Level.WARNING, variant.getKey() + " has no VEP annotations.");
                continue;
            }

            //loop over annotations
            for (VEPAnnotationv75 annotation : vepAnnotations.get(variant.getKey())){

                if (!annotationNodes.containsKey(variant.getKey())){
                    annotationNodes.put(variant.getKey(), new HashMap<String, Node>());
                }

                //add symbols
                if (annotation.getSymbol() != null && !annotation.getSymbol().equals("") && annotation.getSymbolSource().equals("HGNC")){

                    if (!symbolNodes.containsKey(annotation.getSymbol())){
                        symbolNodes.put(annotation.getSymbol(), Neo4j.matchOrCreateUniqueNode(graphDb, Neo4j.getSymbolLabel(), "SymbolId", annotation.getSymbol()));
                    }

                    //link variant and symbol
                    Neo4j.createRelationship(graphDb, variant.getValue(), symbolNodes.get(annotation.getSymbol()), Neo4j.getHasInSymbolRelationship(), properties, false);

                }

                //add feature
                if (annotation.getFeature() != null && !annotation.getFeature().equals("") && !featureNodes.containsKey(annotation.getFeature())){
                    ArrayList<Node> nodes = Neo4j.getNodes(graphDb, Neo4j.getFeatureLabel(), "FeatureId", annotation.getFeature());

                    if (nodes.size() == 0) {

                        if(annotation.getFeature() != null) properties.put("FeatureId", annotation.getFeature());
                        if(annotation.getFeatureType() != null) properties.put("FeatureType", annotation.getFeatureType());
                        if(annotation.getStrand() != null){
                            if (annotation.getStrand().equals("1")){
                                properties.put("Strand", true);
                            } else if (annotation.getStrand().equals("-1")){
                                properties.put("Strand", false);
                            }
                        }
                        if(annotation.getExon() != null) properties.put("TotalExons", Short.parseShort(annotation.getExon().split("/")[1]));

                        featureNodes.put(annotation.getFeature(), Neo4j.addNode(graphDb, Neo4j.getFeatureLabel(), properties));
                        properties.clear();

                        //add canonical label
                        if (annotation.getCanonical() != null && annotation.getCanonical().equals("YES")){
                            Neo4j.addNodeLabel(graphDb, featureNodes.get(annotation.getFeature()), Neo4j.getCanonicalLabel());
                        }

                    } else {
                        featureNodes.put(annotation.getFeature(), nodes.get(0));
                    }

                }

                //add functional annotations
                if(annotation.getHgvsCoding() != null) properties.put("HGVSc", annotation.getHgvsCoding());
                if(annotation.getHgvsProtein() != null) properties.put("HGVSp", annotation.getHgvsProtein());
                if(annotation.getExon() != null) properties.put("Exon", Short.parseShort(annotation.getExon().split("/")[0]));
                if(annotation.getIntron() != null) properties.put("Intron", Short.parseShort(annotation.getIntron().split("/")[0]));
                if(annotation.getSift() != null) properties.put("Sift", annotation.getSift());
                if(annotation.getPolyphen() != null) properties.put("Polyphen", annotation.getPolyphen());

                annotationNodes.get(variant.getKey()).put(annotation.getFeature(), Neo4j.addNode(graphDb, Neo4j.getAnnotationLabel(), properties));
                properties.clear();

                //link consequences
                if (annotation.getConsequences().size() > 0){
                    for (String consequence : annotation.getConsequences()){
                        Neo4j.createRelationship(graphDb, newVariantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), DynamicRelationshipType.withName("HAS_" + consequence.toUpperCase() + "_CONSEQUENCE"), properties, false);
                    }
                } else {
                    Neo4j.createRelationship(graphDb, newVariantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), Neo4j.getHasUnknownConsequenceRelationship(), properties, false);
                }

                //add in feature relationship
                if (annotationNodes.get(variant.getKey()).containsKey(annotation.getFeature()) && featureNodes.containsKey(annotation.getFeature())){
                    Neo4j.createRelationship(graphDb, annotationNodes.get(variant.getKey()).get(annotation.getFeature()), featureNodes.get(annotation.getFeature()), Neo4j.getHasInFeatureRelationship(), properties, false);
                }

                //add in symbol relationship
                if (symbolNodes.containsKey(annotation.getSymbol()) && featureNodes.containsKey(annotation.getFeature())){
                    Neo4j.createRelationship(graphDb, symbolNodes.get(annotation.getSymbol()), featureNodes.get(annotation.getFeature()), DynamicRelationshipType.withName("HAS_" + annotation.getBiotype().toUpperCase() + "_BIOTYPE"), properties, false);
                }

            }
        }

    }
    public void addGenePanels() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Adding virtual panels ...");

        HashMap<String, Object> properties = new HashMap<>();
        Node symbolNode, virtualPanel;

        //add BC
        properties.put("PanelName", "Breast Cancer");
        virtualPanel = Neo4j.addNode(graphDb, Neo4j.getVirtualPanelLabel(), properties);
        properties.clear();

        symbolNode = Neo4j.matchOrCreateUniqueNode(graphDb, Neo4j.getSymbolLabel(), "SymbolId", "BRCA1");
        Neo4j.createRelationship(graphDb, virtualPanel, symbolNode, Neo4j.getHasContainsSymbol(), properties, true);

        symbolNode = Neo4j.matchOrCreateUniqueNode(graphDb, Neo4j.getSymbolLabel(), "SymbolId", "BRCA2");
        Neo4j.createRelationship(graphDb, virtualPanel, symbolNode, Neo4j.getHasContainsSymbol(), properties, true);

        properties.put("Date", System.currentTimeMillis());
        Neo4j.createRelationship(graphDb, virtualPanel, userNode, Neo4j.getHasDesignedBy(), properties, true);
        properties.clear();

        //add TS
        properties.put("PanelName", "Tuberous Sclerosis");
        virtualPanel = Neo4j.addNode(graphDb, Neo4j.getVirtualPanelLabel(), properties);
        properties.clear();

        symbolNode = Neo4j.matchOrCreateUniqueNode(graphDb, Neo4j.getSymbolLabel(), "SymbolId", "TSC1");
        Neo4j.createRelationship(graphDb, virtualPanel, symbolNode, Neo4j.getHasContainsSymbol(), properties, true);

        symbolNode = Neo4j.matchOrCreateUniqueNode(graphDb, Neo4j.getSymbolLabel(), "SymbolId", "TSC2");
        Neo4j.createRelationship(graphDb, virtualPanel, symbolNode, Neo4j.getHasContainsSymbol(), properties, true);

        properties.put("Date", System.currentTimeMillis());
        Neo4j.createRelationship(graphDb, virtualPanel, userNode, Neo4j.getHasDesignedBy(), properties, true);
        properties.clear();

    }

    private void addVariantNodesAndGenotypeRelationshipsHelper(GenomeVariant genomeVariant, short genotypeQuality, Node runInfoNode, RelationshipType relationshipType){

        HashMap<String, Object> properties = new HashMap<>();

        if (!newVariantNodes.containsKey(genomeVariant) && !existingNodes.containsKey(genomeVariant)) {
            ArrayList<Node> nodes = Neo4j.getNodes(graphDb, Neo4j.getVariantLabel(), "VariantId", genomeVariant.getConcatenatedVariant());

            if (nodes.size() == 0){
                properties.put("VariantId", genomeVariant.getConcatenatedVariant());

                //add population freqs
                if (genomeVariantAnnotations.containsKey(genomeVariant)){
                    for (Map.Entry<String, Object> iter : genomeVariantAnnotations.get(genomeVariant).entrySet()){
                        properties.put(iter.getKey(), iter.getValue());
                    }
                } else {
                    log.log(Level.WARNING, genomeVariant.getConcatenatedVariant() + " has no annotated VCF entry.");
                }

                newVariantNodes.put(genomeVariant, Neo4j.addNode(graphDb, Neo4j.getVariantLabel(), properties));
                properties.clear();

                if (genomeVariant.getContig().equals("X")){
                    Neo4j.addNodeLabel(graphDb, newVariantNodes.get(genomeVariant), Neo4j.getXChromosomeLabel());
                } else if (genomeVariant.getContig().equals("Y")){
                    Neo4j.addNodeLabel(graphDb, newVariantNodes.get(genomeVariant), Neo4j.getYChromosomeLabel());
                } else if (Integer.parseInt(genomeVariant.getContig()) > 0 && Integer.parseInt(genomeVariant.getContig()) < 23){
                    Neo4j.addNodeLabel(graphDb, newVariantNodes.get(genomeVariant), Neo4j.getAutoChromosomeLabel());
                }

                //create genotype relationship
                properties.put("Quality", genotypeQuality);
                Neo4j.createRelationship(graphDb, runInfoNode, newVariantNodes.get(genomeVariant), relationshipType, properties, false);
                properties.clear();

            } else {
                existingNodes.put(genomeVariant, nodes.get(0));

                //create genotype relationship
                properties.put("Quality", genotypeQuality);
                Neo4j.createRelationship(graphDb, runInfoNode, existingNodes.get(genomeVariant), relationshipType, properties, false);
                properties.clear();
            }

        }

    }

    public void shutdownDatabase(){
        log.log(Level.INFO, "Shutting down database ...");
        Neo4j.shutdownDatabase(graphDb);
    }

}