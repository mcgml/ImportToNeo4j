package nhs.genetics.cardiff;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.neo4j.graphdb.*;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;

import java.io.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by ml on 23/06/15.
 */

public class VariantDatabase {
    private static final Logger log = Logger.getLogger(VariantDatabase.class.getName());

    private File dbPath;
    private GraphDatabaseService graphDb;
    private VCFFileReader vcfFileReader;
    private HashMap<GenomeVariant, Node> addedVariantNodes = new HashMap<>(); //new variants added during this session
    private HashMap<String, Node> runInfoNodes = new HashMap<>(); //analyses added during this session

    //DB model
    private static Label sampleLabel = Label.label("Sample");
    private static Label variantLabel = Label.label("Variant");
    private static Label autosomeLabel = Label.label("Autosome");
    private static Label xChromLabel = Label.label("X");
    private static Label yChromLabel = Label.label("Y");
    private static Label mtChromLabel = Label.label("MT");
    private static Label snpLabel = Label.label("Snp");
    private static Label indelLabel = Label.label("Indel");
    private static Label annotationLabel = Label.label("Annotation");
    private static Label symbolLabel = Label.label("Symbol");
    private static Label canonicalLabel = Label.label("Canonical");
    private static Label featureLabel = Label.label("Feature");
    private static Label runInfoLabel = Label.label("RunInfo");
    private static Label virtualPanelLabel = Label.label("VirtualPanel");
    private static Label userLabel = Label.label("User");
    private static Label featurePreferenceLabel = Label.label("FeaturePreference");
    private static Label variantPathogenicityLabel = Label.label("VariantPathogenicity");
    private static Label qualityControlLabel = Label.label("QualityControl");
    private static Label disorderLabel = Label.label("Disorder");
    private static RelationshipType hasHetVariantRelationship = RelationshipType.withName("HAS_HET_VARIANT");
    private static RelationshipType hasHomVariantRelationship = RelationshipType.withName("HAS_HOM_VARIANT");
    private static RelationshipType inSymbolRelationship = RelationshipType.withName("IN_SYMBOL");
    private static RelationshipType inFeatureRelationship = RelationshipType.withName("IN_FEATURE");
    private static RelationshipType hasUnknownConsequenceRelationship = RelationshipType.withName("HAS_UNKNOWN_CONSEQUENCE");
    private static RelationshipType hasAnalysisRelationship = RelationshipType.withName("HAS_ANALYSIS");
    private static RelationshipType designedByRelationship = RelationshipType.withName("DESIGNED_BY");
    private static RelationshipType containsSymbolRelationship = RelationshipType.withName("CONTAINS_SYMBOL");
    private static RelationshipType hasProteinCodingBiotypeRelationship = RelationshipType.withName("HAS_PROTEIN_CODING_BIOTYPE");
    private static RelationshipType hasUserEventRelationship = RelationshipType.withName("HAS_USER_EVENT");
    private static RelationshipType addedByRelationship = RelationshipType.withName("ADDED_BY");
    private static RelationshipType authorisedByRelationship = RelationshipType.withName("AUTHORISED_BY");
    private static RelationshipType rejectedByRelationship = RelationshipType.withName("REJECTED_BY");
    private static RelationshipType hasAssociatedSymbol = RelationshipType.withName("HAS_ASSOCIATED_SYMBOL");

    //todo add mutation taster
    //todo add splicing tools

    //population frequencies
    public enum exacPopulation {
        AFR, AMR, EAS, NFE, SAS, FIN, OTH
    }

    public enum kGPhase3Population {
        AFR, AMR, EAS, EUR, SAS
    }

    public VariantDatabase(VCFFileReader vcfFileReader, File dbPath){
        this.vcfFileReader = vcfFileReader;
        this.dbPath = dbPath;
    }

    public void startDatabase(){
        log.log(Level.INFO, "Starting database ...");

        graphDb = new GraphDatabaseFactory().newEmbeddedDatabase(dbPath);
        Neo4j.registerShutdownHook(graphDb);
    }

    //new database
    public void createIndexes() {
        log.log(Level.INFO, "Adding constraints ...");

        //todo optimise
        Neo4j.createConstraint(graphDb, sampleLabel, "sampleId");
        Neo4j.createIndex(graphDb, runInfoLabel, "worklistId");
        Neo4j.createIndex(graphDb, runInfoLabel, "seqId");
        Neo4j.createConstraint(graphDb, runInfoLabel, "analysisId");
        Neo4j.createConstraint(graphDb, variantLabel, "variantId");
        Neo4j.createConstraint(graphDb, featureLabel, "featureId");
        Neo4j.createConstraint(graphDb, symbolLabel, "symbolId");
        Neo4j.createConstraint(graphDb, virtualPanelLabel, "virtualPanelId");
        Neo4j.createConstraint(graphDb, userLabel, "userId");
        Neo4j.createConstraint(graphDb, disorderLabel, "disorder");

    }

    //import genotype VCF
    public void addSampleAndRunInfoNodes() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Adding sample and run info nodes ...");

        HashMap<String, Object> properties = new HashMap<>();
        HashMap<String, String> keyValuePairs = new HashMap<>();
        Set<VCFHeaderLine> metaLines = vcfFileReader.getFileHeader().getMetaDataInInputOrder();

        for (VCFHeaderLine line : metaLines){
            if (line.getKey().equals("SAMPLE")){

                //split out key value pairs
                for (String keyValuePair : line.getValue().split(",")){
                    String[] keyValue = keyValuePair.split("=");
                    keyValuePairs.put(keyValue[0].replace("<", ""), keyValue[1].replace(">", ""));
                }

                //add sample
                Node sampleNode = Neo4j.matchOrCreateUniqueNode(graphDb, sampleLabel, "sampleId", keyValuePairs.get("ID"));

                properties.put("tissue", keyValuePairs.get("Tissue"));
                Neo4j.addNodeProperties(graphDb, sampleNode, properties);
                properties.clear();

                //add run info
                properties.put("worklistId", keyValuePairs.get("WorklistId"));
                properties.put("seqId", keyValuePairs.get("SeqId"));
                properties.put("analysisId", keyValuePairs.get("WorklistId") + "_" + keyValuePairs.get("ID") + "_" + keyValuePairs.get("SeqId"));
                properties.put("assay", keyValuePairs.get("Assay"));
                properties.put("pipelineName", keyValuePairs.get("PipelineName"));
                properties.put("pipelineVersion", Integer.parseInt(keyValuePairs.get("PipelineVersion")));
                properties.put("remoteBamFilePath", keyValuePairs.get("RemoteBamFilePath"));
                properties.put("remoteVcfFilePath", keyValuePairs.get("RemoteVcfFilePath"));

                Node runInfoNode = Neo4j.addNode(graphDb, runInfoLabel, properties);
                properties.clear();

                //link sample and runInfo
                Neo4j.createRelationship(graphDb, sampleNode, runInfoNode, hasAnalysisRelationship, null);
                runInfoNodes.put(keyValuePairs.get("ID"), runInfoNode);

                keyValuePairs.clear();

            }
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
                    if (genotype.isNoCall() || genotype.isHomRef() || genotype.isFiltered()){
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

                        addVariantAndGenotype(genomeVariant, (short) genotype.getGQ(), runInfoNodes.get(genotype.getSampleName()), hasHomVariantRelationship);

                    } else if (genotype.isHet()){

                        genomeVariant = new GenomeVariant(variantContext.getContig(), variantContext.getStart(), variantContext.getReference().getBaseString(), genotype.getAlleles().get(1).getBaseString());
                        genomeVariant.convertToMinimalRepresentation();

                        addVariantAndGenotype(genomeVariant, (short) genotype.getGQ(), runInfoNodes.get(genotype.getSampleName()), hasHetVariantRelationship);

                        if (genotype.isHetNonRef()){

                            genomeVariant = new GenomeVariant(variantContext.getContig(), variantContext.getStart(), variantContext.getReference().getBaseString(), genotype.getAlleles().get(0).getBaseString());
                            genomeVariant.convertToMinimalRepresentation();

                            addVariantAndGenotype(genomeVariant, (short) genotype.getGQ(), runInfoNodes.get(genotype.getSampleName()), hasHetVariantRelationship);
                        }

                    } else {
                        throw new InvalidPropertiesFormatException("Inheritance unknown: " + variantContext.toString());
                    }

                }

            }

        }

    }

    public void writeNewVariantsToVCF(){
        log.log(Level.INFO, "Writing imported variants to VCF.");

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
            properties.put("quality", genotypeQuality);
            Neo4j.createRelationship(graphDb, runInfoNode, addedVariantNodes.get(genomeVariant), relationshipType, properties);

        } catch(IllegalArgumentException | NullPointerException absentNodeException){

            properties.clear();

            try {

                //get variant node
                Node variantNode = Neo4j.getNodes(graphDb, variantLabel, "variantId", genomeVariant.getConcatenatedVariant()).get(0);

                //create genotype relationship
                properties.put("quality", genotypeQuality);
                Neo4j.createRelationship(graphDb, runInfoNode, variantNode, relationshipType, properties);

            } catch (IndexOutOfBoundsException indexOutOfBoundsException) {

                properties.clear();

                //add new variant
                properties.put("variantId", genomeVariant.getConcatenatedVariant());
                Node variantNode = Neo4j.addNode(graphDb, variantLabel, properties);
                properties.clear();

                if (genomeVariant.getContig().equals("X")) {
                    Neo4j.addNodeLabel(graphDb, variantNode, xChromLabel);
                } else if (genomeVariant.getContig().equals("Y")) {
                    Neo4j.addNodeLabel(graphDb, variantNode, yChromLabel);
                } else if (Integer.parseInt(genomeVariant.getContig()) > 0 && Integer.parseInt(genomeVariant.getContig()) < 23) {
                    Neo4j.addNodeLabel(graphDb, variantNode, autosomeLabel);
                }

                if (genomeVariant.isSnp()) Neo4j.addNodeLabel(graphDb, variantNode, snpLabel);
                if (genomeVariant.isIndel()) Neo4j.addNodeLabel(graphDb, variantNode, indelLabel);

                //create genotype relationship
                properties.put("quality", genotypeQuality);
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
            
            Node variantNode = Neo4j.getNodes(graphDb, variantLabel, "variantId", variantId).get(0);

            //add dbSNP Id
            if (variantContext.getID() != null && !variantContext.getID().equals("") && !variantContext.getID().equals(".")){

                properties.put("dbSnpId", variantContext.getID());
                Neo4j.addNodeProperties(graphDb, variantNode, properties);

                properties.clear();
            }

            addVepAnnotations(variantNode, variantContext);
            addPopulationFrequencies(variantNode, variantContext);
            addConservationScores(variantNode, variantContext);

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

            if (!filterVepAnnotation(vepAnnotationv82)){
                vepAnnotations.add(vepAnnotationv82);
            }

        } catch (ClassCastException e) {

            //multiple annotations
            for (String annotation : (ArrayList<String>) variantContext.getAttribute("CSQ")) {

                VEPAnnotationv82 vepAnnotationv82 = new VEPAnnotationv82(annotation);
                vepAnnotationv82.parseAnnotation();

                if (!filterVepAnnotation(vepAnnotationv82)) {
                    vepAnnotations.add(vepAnnotationv82);
                }

            }

        }

        //loop over annotations
        for (VEPAnnotationv82 annotation : vepAnnotations) {

            symbolNode = null;
            featureNode = null;
            annotationNode = null;

            //add symbol
            if (annotation.getSymbol() != null && !annotation.getSymbol().equals("")) {
                symbolNode = Neo4j.matchOrCreateUniqueNode(graphDb, symbolLabel, "symbolId", annotation.getSymbol()); //add symbol
                Neo4j.createRelationship(graphDb, variantNode, symbolNode, inSymbolRelationship, properties); //link variant and symbol
            }

            //add feature
            if (annotation.getFeature() != null && !annotation.getFeature().equals("")) {
                featureNode = Neo4j.matchOrCreateUniqueNode(graphDb, featureLabel, "featureId", annotation.getFeature()); //add feature

                if (annotation.getFeature() != null) properties.put("featureId", annotation.getFeature());
                if (annotation.getFeatureType() != null) properties.put("featureType", annotation.getFeatureType());
                if (annotation.getCcds() != null) properties.put("ccdsId", annotation.getCcds());
                if (annotation.getStrand() == 1) {
                    properties.put("strand", true);
                } else if (annotation.getStrand() == -1) {
                    properties.put("strand", false);
                }
                if (annotation.getExon() != null) properties.put("totalExons", Short.parseShort(annotation.getExon().split("/")[1]));

                Neo4j.addNodeProperties(graphDb, featureNode, properties);
                properties.clear();

                if (annotation.isCanonical()) {
                    Neo4j.addNodeLabel(graphDb, featureNode, canonicalLabel);
                }
            }

            //add annotation
            if (annotation.getHgvsc() != null) properties.put("hgvsc", annotation.getHgvsc());
            if (annotation.getHgvsp() != null) properties.put("hgvsp", annotation.getHgvsp());
            if (annotation.getExon() != null) properties.put("exon", annotation.getExon().split("/")[0]); //must remain as string, can be given as range i.e. 1-2
            if (annotation.getIntron() != null) properties.put("intron", annotation.getIntron().split("/")[0]); //must remain as string, can be given as range i.e. 1-2
            if (annotation.getSift() != null) properties.put("sift", annotation.getSift());
            if (annotation.getPolyPhen() != null) properties.put("polyphen", annotation.getPolyPhen());
            if (annotation.getCodons() != null) properties.put("codons", annotation.getCodons());

            //add protein domains Pfam_domain
            if (annotation.getDomains().containsKey("Pfam_domain")){
                properties.put("pfamDomain", annotation.getDomains().get("Pfam_domain").toArray(new String[annotation.getDomains().get("Pfam_domain").size()]));
            }

            //add protein domains hmmpanther
            if (annotation.getDomains().containsKey("hmmpanther")){
                properties.put("hmmPanther", annotation.getDomains().get("hmmpanther").toArray(new String[annotation.getDomains().get("hmmpanther").size()]));
            }

            //add protein domains PROSITE_profiles && PROSITE_patterns
            if (annotation.getDomains().containsKey("PROSITE_profiles") || annotation.getDomains().containsKey("PROSITE_patterns")){

                //combine Prosite numbers
                HashSet<String> temp = new HashSet<>();
                if (annotation.getDomains().containsKey("PROSITE_profiles")) temp.addAll(annotation.getDomains().get("PROSITE_profiles"));
                if (annotation.getDomains().containsKey("PROSITE_patterns")) temp.addAll(annotation.getDomains().get("PROSITE_patterns"));

                properties.put("prosite", temp.toArray(new String[temp.size()]));
            }

            //add protein domains Superfamily_domains
            if (annotation.getDomains().containsKey("Superfamily_domains")){
                properties.put("superfamilyDomains", annotation.getDomains().get("Superfamily_domains").toArray(new String[annotation.getDomains().get("Superfamily_domains").size()]));
            }

            annotationNode = Neo4j.addNode(graphDb, annotationLabel, properties);
            properties.clear();

            //link consequences
            if (annotation.getConsequences().size() > 0) {
                for (String consequence : annotation.getConsequences()) {
                    Neo4j.createRelationship(graphDb, variantNode, annotationNode, RelationshipType.withName("HAS_" + consequence.toUpperCase() + "_CONSEQUENCE"), properties);
                }
            } else {
                Neo4j.createRelationship(graphDb, variantNode, annotationNode, hasUnknownConsequenceRelationship, properties);
            }

            //add in feature relationship
            if (annotationNode != null && featureNode != null) {
                Neo4j.createRelationship(graphDb, annotationNode, featureNode, inFeatureRelationship, properties);
            }

            //add in symbol relationship
            if (symbolNode != null && featureNode != null) {
                Neo4j.createRelationship(graphDb, symbolNode, featureNode, RelationshipType.withName("HAS_" + annotation.getBiotype().toUpperCase() + "_BIOTYPE"), properties);
            }

        }

    }

    private void addPopulationFrequencies(Node variantNode, VariantContext variantContext){

        int minimumAllelesForAFCalculation = 120;
        HashMap<String, Object> properties = new HashMap<>();

        // 1000 genomes phase 3
        for (kGPhase3Population populationFrequency : kGPhase3Population.values()){
            if (variantContext.getAttribute("kGPhase3." + populationFrequency.toString() + "_AF") != null && !variantContext.getAttribute("kGPhase3." + populationFrequency.toString() + "_AF").equals(".")) {
                properties.put("kGPhase3" + populationFrequency.toString() + "Af", Float.parseFloat((String) variantContext.getAttribute("kGPhase3." + populationFrequency.toString() + "_AF")));
            }
        }

        // Exome aggregation consortium
        for (exacPopulation populationFrequency : exacPopulation.values()){
            if (variantContext.getAttribute("exac.AC_" + populationFrequency.toString()) != null && !variantContext.getAttribute("exac.AC_" + populationFrequency.toString()).equals(".")
                    && variantContext.getAttribute("exac.AN_" + populationFrequency.toString()) != null && !variantContext.getAttribute("exac.AN_" + populationFrequency.toString()).equals(".")
                    && Integer.parseInt((String) variantContext.getAttribute("exac.AN_" + populationFrequency.toString())) > minimumAllelesForAFCalculation) {
                properties.put("exac" + populationFrequency.toString() + "Af", Float.parseFloat((String) variantContext.getAttribute("exac.AC_" + populationFrequency.toString())) / Float.parseFloat((String) variantContext.getAttribute("exac.AN_" + populationFrequency.toString())));
            }
        }

        Neo4j.addNodeProperties(graphDb, variantNode, properties);

    }

    private void addConservationScores(Node variantNode, VariantContext variantContext){
        HashMap<String, Object> properties = new HashMap<>();

        if (variantContext.getAttribute("GERP") != null && !variantContext.getAttribute("GERP").equals(".")) {
            properties.put("gerp", Float.parseFloat((String) variantContext.getAttribute("GERP")));
        }
        if (variantContext.getAttribute("phastCons") != null && !variantContext.getAttribute("phastCons").equals(".")) {
            properties.put("phastCons", Float.parseFloat((String) variantContext.getAttribute("phastCons")));
        }
        if (variantContext.getAttribute("phyloP") != null && !variantContext.getAttribute("phyloP").equals(".")) {
            properties.put("phyloP", Float.parseFloat((String) variantContext.getAttribute("phyloP")));
        }

        Neo4j.addNodeProperties(graphDb, variantNode, properties);

    }

    private static boolean filterVepAnnotation(VEPAnnotationv82 vepAnnotationv82) {

        //check biotype
        if (vepAnnotationv82.getBiotype() == null || !vepAnnotationv82.getBiotype().equals("protein_coding")) {
            return true;
        }

        //check symbol source
        if (vepAnnotationv82.getSymbolSource() == null || !vepAnnotationv82.getSymbolSource().equals("HGNC")) {
            return true;
        }

        return false;
    }

    public void shutdownDatabase(){
        log.log(Level.INFO, "Shutting down database ...");
        Neo4j.shutdownDatabase(graphDb);
    }

    public static Label getSampleLabel() {
        return sampleLabel;
    }

    public static Label getVariantLabel() {
        return variantLabel;
    }

    public static Label getAutosomeLabel() {
        return autosomeLabel;
    }

    public static Label getxChromLabel() {
        return xChromLabel;
    }

    public static Label getyChromLabel() {
        return yChromLabel;
    }

    public static Label getMtChromLabel() {
        return mtChromLabel;
    }

    public static Label getSnpLabel() {
        return snpLabel;
    }

    public static Label getIndelLabel() {
        return indelLabel;
    }

    public static Label getAnnotationLabel() {
        return annotationLabel;
    }

    public static Label getSymbolLabel() {
        return symbolLabel;
    }

    public static Label getDisorderLabel() {
        return disorderLabel;
    }

    public static Label getCanonicalLabel() {
        return canonicalLabel;
    }

    public static Label getFeatureLabel() {
        return featureLabel;
    }

    public static Label getRunInfoLabel() {
        return runInfoLabel;
    }

    public static Label getVirtualPanelLabel() {
        return virtualPanelLabel;
    }

    public static Label getUserLabel() {
        return userLabel;
    }

    public static Label getFeaturePreferenceLabel() {
        return featurePreferenceLabel;
    }

    public static Label getVariantPathogenicityLabel() {
        return variantPathogenicityLabel;
    }

    public static RelationshipType getHasHetVariantRelationship() {
        return hasHetVariantRelationship;
    }

    public static RelationshipType getHasHomVariantRelationship() {
        return hasHomVariantRelationship;
    }

    public static RelationshipType getInSymbolRelationship() {
        return inSymbolRelationship;
    }

    public static RelationshipType getInFeatureRelationship() {
        return inFeatureRelationship;
    }

    public static RelationshipType getHasUnknownConsequenceRelationship() {
        return hasUnknownConsequenceRelationship;
    }

    public static RelationshipType getHasAnalysisRelationship() {
        return hasAnalysisRelationship;
    }

    public static RelationshipType getDesignedByRelationship() {
        return designedByRelationship;
    }

    public static RelationshipType getContainsSymbolRelationship() {
        return containsSymbolRelationship;
    }

    public static RelationshipType getHasProteinCodingBiotypeRelationship() {
        return hasProteinCodingBiotypeRelationship;
    }

    public static RelationshipType getHasUserEventRelationship() {
        return hasUserEventRelationship;
    }

    public static RelationshipType getAddedByRelationship() {
        return addedByRelationship;
    }

    public static RelationshipType getAuthorisedByRelationship() {
        return authorisedByRelationship;
    }

    public static RelationshipType getRejectedByRelationship() {
        return rejectedByRelationship;
    }

    public static RelationshipType getHasAssociatedSymbol() {
        return hasAssociatedSymbol;
    }

    public static Label getQualityControlLabel() {
        return qualityControlLabel;
    }
}