package nhs.genetics.cardiff;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.neo4j.graphdb.*;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;

import java.io.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by ml on 23/06/15.
 */

//TODO retrieve pubmed abstracts (web ui)
//TODO check alamut for extra functionality
//todo add ddd afs
//todo add clinvar

public class VariantDatabase {
    private static final Logger log = Logger.getLogger(VariantDatabase.class.getName());

    private int pipelineVersion, minimumAlellesForAFCalculation = 120;
    private String worklistId, runId, supplierPanelName, pipelineName;
    private File dbPath;
    private GraphDatabaseService graphDb;
    private VCFFileReader vcfFileReader;
    private HashMap<GenomeVariant, Node> addedVariantNodes = new HashMap<>(); //new variants added during this session
    private HashMap<String, Node> runInfoNodes = new HashMap<>(); //analyses added during this session
    private Node userNode;

    public VariantDatabase(VCFFileReader vcfFileReader, File dbPath){
        this.vcfFileReader = vcfFileReader;
        this.dbPath = dbPath;
    }
    public void startDatabase(){
        log.log(Level.INFO, "Starting database ...");

        graphDb = new GraphDatabaseFactory()
                .newEmbeddedDatabaseBuilder(dbPath)
                .newGraphDatabase();

        Neo4j.registerShutdownHook(graphDb);
    }

    //new database
    public void createIndexes() {
        log.log(Level.INFO, "Adding constraints ...");

        //todo check these
        Neo4j.createConstraint(graphDb, Neo4j.getSampleLabel(), "SampleId");
        Neo4j.createIndex(graphDb, Neo4j.getRunInfoLabel(), "WorklistId");
        Neo4j.createIndex(graphDb, Neo4j.getRunInfoLabel(), "RunId");
        Neo4j.createConstraint(graphDb, Neo4j.getRunInfoLabel(), "AnalysisId");
        Neo4j.createConstraint(graphDb, Neo4j.getVariantLabel(), "VariantId");
        Neo4j.createConstraint(graphDb, Neo4j.getMtChromosomeLabel(), "VariantId");
        Neo4j.createConstraint(graphDb, Neo4j.getXChromosomeLabel(), "VariantId");
        Neo4j.createConstraint(graphDb, Neo4j.getYChromosomeLabel(), "VariantId");
        Neo4j.createConstraint(graphDb, Neo4j.getAutoChromosomeLabel(), "VariantId");
        Neo4j.createIndex(graphDb, Neo4j.getVariantLabel(), "dbSNPId");
        Neo4j.createConstraint(graphDb, Neo4j.getFeatureLabel(), "FeatureId");
        Neo4j.createConstraint(graphDb, Neo4j.getCanonicalLabel(), "FeatureId");
        Neo4j.createConstraint(graphDb, Neo4j.getSymbolLabel(), "SymbolId");
        Neo4j.createConstraint(graphDb, Neo4j.getVirtualPanelLabel(), "VirtualPanelId");
        Neo4j.createConstraint(graphDb, Neo4j.getUserLabel(), "UserId");

    }
    public void addUsers() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Adding users ...");
        userNode = Neo4j.matchOrCreateUniqueNode(graphDb, Neo4j.getUserLabel(), "UserId", "ml");
    }
    public void addVirtualPanels() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Adding virtual panels ...");

        HashMap<String, Object> properties = new HashMap<>();
        Node symbolNode, virtualPanel;

        //add BC
        properties.put("VirtualPanelName", "Breast Cancer");
        virtualPanel = Neo4j.addNode(graphDb, Neo4j.getVirtualPanelLabel(), properties);
        properties.clear();

        symbolNode = Neo4j.matchOrCreateUniqueNode(graphDb, Neo4j.getSymbolLabel(), "SymbolId", "BRCA1");
        Neo4j.createRelationship(graphDb, virtualPanel, symbolNode, Neo4j.getHasContainsSymbol(), properties);

        symbolNode = Neo4j.matchOrCreateUniqueNode(graphDb, Neo4j.getSymbolLabel(), "SymbolId", "BRCA2");
        Neo4j.createRelationship(graphDb, virtualPanel, symbolNode, Neo4j.getHasContainsSymbol(), properties);

        properties.put("Date", System.currentTimeMillis());
        Neo4j.createRelationship(graphDb, virtualPanel, userNode, Neo4j.getHasDesignedBy(), properties);
        properties.clear();

        //add TS
        properties.put("VirtualPanelName", "Tuberous Sclerosis");
        virtualPanel = Neo4j.addNode(graphDb, Neo4j.getVirtualPanelLabel(), properties);
        properties.clear();

        symbolNode = Neo4j.matchOrCreateUniqueNode(graphDb, Neo4j.getSymbolLabel(), "SymbolId", "TSC1");
        Neo4j.createRelationship(graphDb, virtualPanel, symbolNode, Neo4j.getHasContainsSymbol(), properties);

        symbolNode = Neo4j.matchOrCreateUniqueNode(graphDb, Neo4j.getSymbolLabel(), "SymbolId", "TSC2");
        Neo4j.createRelationship(graphDb, virtualPanel, symbolNode, Neo4j.getHasContainsSymbol(), properties);

        properties.put("Date", System.currentTimeMillis());
        Neo4j.createRelationship(graphDb, virtualPanel, userNode, Neo4j.getHasDesignedBy(), properties);
        properties.clear();

    }

    //import genotype VCF
    public void populateRunMetaData(){
        log.log(Level.INFO, "Extracting metadata ...");

        worklistId = "K15-0000";
        runId = "AFWJ3";
        supplierPanelName = "Illumina TruSight Cancer";
        pipelineName = "IlluminaTruSightCancer";
        pipelineVersion = 4;

    }
    public void addSampleAndRunInfoNodes() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Adding sample and run info nodes ...");

        HashMap<String, Object> properties = new HashMap<>();
        ArrayList<String> sampleIds = vcfFileReader.getFileHeader().getSampleNamesInOrder();

        for (short n = 0; n < sampleIds.size(); ++n){

            //add sample
            Node sampleNode = Neo4j.matchOrCreateUniqueNode(graphDb, Neo4j.getSampleLabel(), "SampleId", sampleIds.get(n));

            properties.put("SampleType", "Peripheral Blood"); //todo
            Neo4j.addNodeProperties(graphDb, sampleNode, properties);
            properties.clear();

            //add run info
            properties.put("WorklistId", worklistId);
            properties.put("RunId", runId);
            properties.put("SampleNo", n);
            properties.put("AnalysisId", worklistId + "_" + n + "_" + runId);
            properties.put("SupplierPanelName", supplierPanelName);
            properties.put("PipelineName", pipelineName);
            properties.put("PipelineVersion", pipelineVersion);
            properties.put("RemoteBamFilePath", "bam.bam");//todo
            properties.put("RemoteVcfFilePath", "vcf.vcf");//todo

            Node runInfoNode = Neo4j.addNode(graphDb, Neo4j.getRunInfoLabel(), properties);
            properties.clear();

            //link sample and runInfo
            Neo4j.createRelationship(graphDb, sampleNode, runInfoNode, Neo4j.getHasAnalysisRelationship(), properties);
            runInfoNodes.put(sampleIds.get(n), runInfoNode);
        }

    }
    public void importVariants() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Importing variants ...");

        GenomeVariant genomeVariant;
        Iterator<VariantContext> variantContextIterator = vcfFileReader.iterator();

        //read variant VCF file
        while (variantContextIterator.hasNext()) {
            VariantContext variantContext = variantContextIterator.next();

            //skip filtered and non-variant loci
            if (!variantContext.isFiltered() && variantContext.isVariant()){
                Iterator<Genotype> genotypeIterator = variantContext.getGenotypes().iterator();

                //read genotypes
                while (genotypeIterator.hasNext()) {
                    Genotype genotype = genotypeIterator.next();

                    //skip no-calls, hom-refs,  mixed genotypes or alleles covered by nearby indels
                    if (genotype.isNoCall() || genotype.isHomRef()){
                        continue;
                    }
                    if (genotype.isMixed()){
                        log.log(Level.WARNING, genotype.getSampleName() + ": " + variantContext.getContig() + " " + variantContext.getStart() + " " + variantContext.getReference() + variantContext.getAlternateAlleles().toString() + " has mixed genotype ( " + genotype.getGenotypeString() + " ) and could not be added.");
                        continue;
                    }
                    if (genotype.getPloidy() != 2 || genotype.getAlleles().size() != 2) {
                        throw new InvalidPropertiesFormatException("Allele " + genotype.getAlleles().toString() + " is not diploid");
                    }
                    if (genotype.getAlleles().get(0).getBaseString().equals("*") || genotype.getAlleles().get(1).getBaseString().equals("*")) {
                        continue;
                    }

                    //add new variants to the DB
                    if (genotype.isHom()){

                        genomeVariant = new GenomeVariant(variantContext.getContig(), variantContext.getStart(), variantContext.getReference().getBaseString(), genotype.getAlleles().get(1).getBaseString());
                        genomeVariant.convertToMinimalRepresentation();

                        addVariantAndGenotype(genomeVariant, (short) genotype.getGQ(), runInfoNodes.get(genotype.getSampleName()), Neo4j.getHasHomVariantRelationship());

                    } else if (genotype.isHet()){

                        genomeVariant = new GenomeVariant(variantContext.getContig(), variantContext.getStart(), variantContext.getReference().getBaseString(), genotype.getAlleles().get(1).getBaseString());
                        genomeVariant.convertToMinimalRepresentation();

                        addVariantAndGenotype(genomeVariant, (short) genotype.getGQ(), runInfoNodes.get(genotype.getSampleName()), Neo4j.getHasHetVariantRelationship());

                        if (genotype.isHetNonRef()){

                            genomeVariant = new GenomeVariant(variantContext.getContig(), variantContext.getStart(), variantContext.getReference().getBaseString(), genotype.getAlleles().get(0).getBaseString());
                            genomeVariant.convertToMinimalRepresentation();

                            addVariantAndGenotype(genomeVariant, (short) genotype.getGQ(), runInfoNodes.get(genotype.getSampleName()), Neo4j.getHasHetVariantRelationship());
                        }

                    } else {
                        throw new InvalidPropertiesFormatException("Inheritance unknown: " + variantContext.toString());
                    }

                }

            }

        }

    }
    public void writeNewVariantsToVCF(){

        try (PrintWriter printWriter = new PrintWriter(new File("imported.vcf"))){

            printWriter.println("##fileformat=VCFv4.1");
            printWriter.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");

            //write out variants
            for (Map.Entry<GenomeVariant, Node> iter : addedVariantNodes.entrySet()){

                printWriter.println
                        (
                                iter.getKey().getContig() + "\t" +
                                        iter.getKey().getPos() + "\t" +
                                        "." + "\t" +
                                        iter.getKey().getRef() + "\t" +
                                        iter.getKey().getAlt() + "\t" +
                                        "." + "\t" +
                                        "." + "\t" +
                                        "."
                        );

            }
        } catch (IOException e){
            log.log(Level.SEVERE, "Could not output variants.");
        }

    }
    private void addVariantAndGenotype(GenomeVariant genomeVariant, short genotypeQuality, Node runInfoNode, RelationshipType relationshipType){
        HashMap<String, Object> properties = new HashMap<>();

        try {

            //create genotype relationship
            properties.put("Quality", genotypeQuality);
            Neo4j.createRelationship(graphDb, runInfoNode, addedVariantNodes.get(genomeVariant), relationshipType, properties);

        } catch(IllegalArgumentException | NullPointerException absentNodeException){

            properties.clear();

            try {

                //get variant node
                Node variantNode = Neo4j.getNodes(graphDb, Neo4j.getVariantLabel(), "VariantId", genomeVariant.getConcatenatedVariant()).get(0);

                //create genotype relationship
                properties.put("Quality", genotypeQuality);
                Neo4j.createRelationship(graphDb, runInfoNode, variantNode, relationshipType, properties);

            } catch (IndexOutOfBoundsException indexOutOfBoundsException) {

                properties.clear();

                //variant not present in graph

                //add new variant
                properties.put("VariantId", genomeVariant.getConcatenatedVariant());
                Node variantNode = Neo4j.addNode(graphDb, Neo4j.getVariantLabel(), properties);
                properties.clear();

                if (genomeVariant.getContig().equals("X")) {
                    Neo4j.addNodeLabel(graphDb, variantNode, Neo4j.getXChromosomeLabel());
                } else if (genomeVariant.getContig().equals("Y")) {
                    Neo4j.addNodeLabel(graphDb, variantNode, Neo4j.getYChromosomeLabel());
                } else if (Integer.parseInt(genomeVariant.getContig()) > 0 && Integer.parseInt(genomeVariant.getContig()) < 23) {
                    Neo4j.addNodeLabel(graphDb, variantNode, Neo4j.getAutoChromosomeLabel());
                }

                //create genotype relationship
                properties.put("Quality", genotypeQuality);
                Neo4j.createRelationship(graphDb, runInfoNode, variantNode, relationshipType, properties);
                properties.clear();

                addedVariantNodes.put(genomeVariant, variantNode);
            }

        }

    }

    //import annotation VCF
    public void importAnnotations() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Importing annotations ...");

        HashMap<String, Object> properties = new HashMap<>();

        Iterator<VariantContext> variantContextIterator = vcfFileReader.iterator();

        //read annotation VCF file
        while (variantContextIterator.hasNext()) {
            VariantContext variantContext = variantContextIterator.next();

            //loop up variant Node
            String variantId = variantContext.getContig() + ":" +
                    variantContext.getStart() +
                    variantContext.getAlleles().get(0).getBaseString() + ">" +
                    variantContext.getAlleles().get(1).getBaseString();

            log.log(Level.INFO, "Adding annotation for: " + variantId);
            Node variantNode = Neo4j.getNodes(graphDb, Neo4j.getVariantLabel(), "VariantId", variantId).get(0);

            //add dbSNP Id
            if (variantContext.getID() != null && !variantContext.getID().equals("") && !variantContext.getID().equals(".")){

                properties.put("dbSNPId", variantContext.getID());
                Neo4j.addNodeProperties(graphDb, variantNode, properties);

                properties.clear();
            }

            addVepAnnotations(variantNode, variantContext);
            addPopulationFrequencies(variantNode, variantContext);

        }
    }
    private void addVepAnnotations(Node variantNode, VariantContext variantContext) throws InvalidPropertiesFormatException {

        HashMap<String, Object> properties = new HashMap<>();
        HashSet<VEPAnnotationv82> vepAnnotations = new HashSet<>();
        Node symbolNode, featureNode, annotationNode;

        //split annotations and make unique
        try {

            //one annotation
            VEPAnnotationv82 vepAnnotationv82 = new VEPAnnotationv82((String) variantContext.getAttribute("CSQ"));
            vepAnnotationv82.parseAnnotation();

            vepAnnotations.add(vepAnnotationv82);

        } catch (ClassCastException e) {

            //multiple annotations
            for (String annotation : (ArrayList<String>) variantContext.getAttribute("CSQ")) {

                VEPAnnotationv82 vepAnnotationv82 = new VEPAnnotationv82(annotation);
                vepAnnotationv82.parseAnnotation();

                vepAnnotations.add(vepAnnotationv82);
            }
        }

        //loop over annotations
        for (VEPAnnotationv82 annotation : vepAnnotations) {

            symbolNode = null;
            featureNode = null;
            annotationNode = null;

            //add symbol
            if (annotation.getSymbol() != null && !annotation.getSymbol().equals("")) {
                symbolNode = Neo4j.matchOrCreateUniqueNode(graphDb, Neo4j.getSymbolLabel(), "SymbolId", annotation.getSymbol()); //add symbol

                if (annotation.getGene() != null && !annotation.getGene().equals("")){
                    properties.put("GeneId", annotation.getGene());
                    Neo4j.addNodeProperties(graphDb, symbolNode, properties);
                    properties.clear();
                }

                Neo4j.createRelationship(graphDb, variantNode, symbolNode, Neo4j.getHasInSymbolRelationship(), properties); //link variant and symbol

            }

            //add feature
            if (annotation.getFeature() != null && !annotation.getFeature().equals("")) {
                featureNode = Neo4j.matchOrCreateUniqueNode(graphDb, Neo4j.getFeatureLabel(), "FeatureId", annotation.getFeature()); //add feature

                if (annotation.getFeature() != null) properties.put("FeatureId", annotation.getFeature());
                if (annotation.getFeatureType() != null) properties.put("FeatureType", annotation.getFeatureType());
                if (annotation.getCcds() != null) properties.put("CCDSId", annotation.getCcds());
                if (annotation.getStrand() == 1) {
                    properties.put("Strand", true);
                } else if (annotation.getStrand() == -1) {
                    properties.put("Strand", false);
                }
                if (annotation.getExon() != null) properties.put("TotalExons", Short.parseShort(annotation.getExon().split("/")[1]));

                Neo4j.addNodeProperties(graphDb, featureNode, properties);
                properties.clear();

                if (annotation.isCanonical()) {
                    Neo4j.addNodeLabel(graphDb, featureNode, Neo4j.getCanonicalLabel());
                }
            }

            //add annotation
            if (annotation.getHgvsCoding() != null) properties.put("HGVSc", annotation.getHgvsCoding());
            if (annotation.getHgvsProtein() != null) properties.put("HGVSp", annotation.getHgvsProtein());
            if (annotation.getExon() != null)
                properties.put("Exon", Short.parseShort(annotation.getExon().split("/")[0]));
            if (annotation.getIntron() != null)
                properties.put("Intron", Short.parseShort(annotation.getIntron().split("/")[0]));
            if (annotation.getSift() != null) properties.put("Sift", annotation.getSift());
            if (annotation.getPolyPhen() != null) properties.put("Polyphen", annotation.getPolyPhen());
            if (annotation.getCodons() != null) properties.put("Codons", annotation.getCodons());

            annotationNode = Neo4j.addNode(graphDb, Neo4j.getAnnotationLabel(), properties);
            properties.clear();

            //link consequences
            if (annotation.getConsequences().size() > 0) {
                for (String consequence : annotation.getConsequences()) {
                    Neo4j.createRelationship(graphDb, variantNode, annotationNode, DynamicRelationshipType.withName("HAS_" + consequence.toUpperCase() + "_CONSEQUENCE"), properties);
                }
            } else {
                Neo4j.createRelationship(graphDb, variantNode, annotationNode, Neo4j.getHasUnknownConsequenceRelationship(), properties);
            }

            //add in feature relationship
            if (annotationNode != null && featureNode != null) {
                Neo4j.createRelationship(graphDb, annotationNode, featureNode, Neo4j.getHasInFeatureRelationship(), properties);
            }

            //add in symbol relationship
            if (symbolNode != null && featureNode != null) {
                Neo4j.createRelationship(graphDb, symbolNode, featureNode, DynamicRelationshipType.withName("HAS_" + annotation.getBiotype().toUpperCase() + "_BIOTYPE"), properties);
            }

        }

    }
    private void addPopulationFrequencies(Node variantNode, VariantContext variantContext){

        HashMap<String, Object> properties = new HashMap<>();

        if (variantContext.getAttribute("1kgp3.EAS_AF") != null && !variantContext.getAttribute("1kgp3.EAS_AF").equals(".")) {
            properties.put("onekGPhase3_EAS_AF", Float.parseFloat((String) variantContext.getAttribute("1kgp3.EAS_AF")));
        }
        if (variantContext.getAttribute("1kgp3.EUR_AF") != null && !variantContext.getAttribute("1kgp3.EUR_AF").equals(".")) {
            properties.put("onekGPhase3_EUR_AF", Float.parseFloat((String) variantContext.getAttribute("1kgp3.EUR_AF")));
        }
        if (variantContext.getAttribute("1kgp3.AFR_AF") != null && !variantContext.getAttribute("1kgp3.AFR_AF").equals(".")) {
            properties.put("onekGPhase3_AFR_AF", Float.parseFloat((String) variantContext.getAttribute("1kgp3.AFR_AF")));
        }
        if (variantContext.getAttribute("1kgp3.AMR_AF") != null && !variantContext.getAttribute("1kgp3.AMR_AF").equals(".")) {
            properties.put("onekGPhase3_AMR_AF", Float.parseFloat((String) variantContext.getAttribute("1kgp3.AMR_AF")));
        }
        if (variantContext.getAttribute("1kgp3.SAS_AF") != null && !variantContext.getAttribute("1kgp3.SAS_AF").equals(".")) {
            properties.put("onekGPhase3_SAS_AF", Float.parseFloat((String) variantContext.getAttribute("1kgp3.SAS_AF")));
        }
        if (variantContext.getAttribute("ExAC.AC_AFR") != null && !variantContext.getAttribute("ExAC.AC_AFR").equals(".") && variantContext.getAttribute("ExAC.AN_AFR") != null && !variantContext.getAttribute("ExAC.AN_AFR").equals(".") && Integer.parseInt((String) variantContext.getAttribute("ExAC.AN_AFR")) > minimumAlellesForAFCalculation) {
            properties.put("ExAC_AFR_AF", Float.parseFloat((String) variantContext.getAttribute("ExAC.AC_AFR")) / Float.parseFloat((String) variantContext.getAttribute("ExAC.AN_AFR")));
        }
        if (variantContext.getAttribute("ExAC.AC_AMR") != null && !variantContext.getAttribute("ExAC.AC_AMR").equals(".") && variantContext.getAttribute("ExAC.AN_AMR") != null && !variantContext.getAttribute("ExAC.AN_AMR").equals(".") && Integer.parseInt((String) variantContext.getAttribute("ExAC.AN_AMR")) > minimumAlellesForAFCalculation) {
            properties.put("ExAC_AMR_AF", Float.parseFloat((String) variantContext.getAttribute("ExAC.AC_AMR")) / Float.parseFloat((String) variantContext.getAttribute("ExAC.AN_AMR")));
        }
        if (variantContext.getAttribute("ExAC.AC_EAS") != null && !variantContext.getAttribute("ExAC.AC_EAS").equals(".") && variantContext.getAttribute("ExAC.AN_EAS") != null && !variantContext.getAttribute("ExAC.AN_EAS").equals(".") && Integer.parseInt((String) variantContext.getAttribute("ExAC.AN_EAS")) > minimumAlellesForAFCalculation) {
            properties.put("ExAC_EAS_AF", Float.parseFloat((String) variantContext.getAttribute("ExAC.AC_EAS")) / Float.parseFloat((String) variantContext.getAttribute("ExAC.AN_EAS")));
        }
        if (variantContext.getAttribute("ExAC.AC_FIN") != null && !variantContext.getAttribute("ExAC.AC_FIN").equals(".") && variantContext.getAttribute("ExAC.AN_FIN") != null && !variantContext.getAttribute("ExAC.AN_FIN").equals(".") && Integer.parseInt((String) variantContext.getAttribute("ExAC.AN_FIN")) > minimumAlellesForAFCalculation) {
            properties.put("ExAC_FIN_AF", Float.parseFloat((String) variantContext.getAttribute("ExAC.AC_FIN")) / Float.parseFloat((String) variantContext.getAttribute("ExAC.AN_FIN")));
        }
        if (variantContext.getAttribute("ExAC.AC_NFE") != null && !variantContext.getAttribute("ExAC.AC_NFE").equals(".") && variantContext.getAttribute("ExAC.AN_NFE") != null && !variantContext.getAttribute("ExAC.AN_NFE").equals(".") && Integer.parseInt((String) variantContext.getAttribute("ExAC.AN_NFE")) > minimumAlellesForAFCalculation) {
            properties.put("ExAC_NFE_AF", Float.parseFloat((String) variantContext.getAttribute("ExAC.AC_NFE")) / Float.parseFloat((String) variantContext.getAttribute("ExAC.AN_NFE")));
        }
        if (variantContext.getAttribute("ExAC.AC_OTH") != null && !variantContext.getAttribute("ExAC.AC_OTH").equals(".") && variantContext.getAttribute("ExAC.AN_OTH") != null && !variantContext.getAttribute("ExAC.AN_OTH").equals(".") && Integer.parseInt((String) variantContext.getAttribute("ExAC.AN_OTH")) > minimumAlellesForAFCalculation) {
            properties.put("ExAC_OTH_AF", Float.parseFloat((String) variantContext.getAttribute("ExAC.AC_OTH")) / Float.parseFloat((String) variantContext.getAttribute("ExAC.AN_OTH")));
        }
        if (variantContext.getAttribute("ExAC.AC_SAS") != null && !variantContext.getAttribute("ExAC.AC_SAS").equals(".") && variantContext.getAttribute("ExAC.AN_SAS") != null && !variantContext.getAttribute("ExAC.AN_SAS").equals(".") && Integer.parseInt((String) variantContext.getAttribute("ExAC.AN_SAS")) > minimumAlellesForAFCalculation) {
            properties.put("ExAC_SAS_AF", Float.parseFloat((String) variantContext.getAttribute("ExAC.AC_SAS")) / Float.parseFloat((String) variantContext.getAttribute("ExAC.AN_SAS")));
        }

        Neo4j.addNodeProperties(graphDb, variantNode, properties);

    }

    public void shutdownDatabase(){
        log.log(Level.INFO, "Shutting down database ...");
        Neo4j.shutdownDatabase(graphDb);
    }

}