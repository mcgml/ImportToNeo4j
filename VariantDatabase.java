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
public class VariantDatabase {

    private static final Logger log = Logger.getLogger(VariantDatabase.class.getName());

    private Label patientLabel = DynamicLabel.label("Patient");
    private Label sampleLabel = DynamicLabel.label("Sample");
    private Label variantLabel = DynamicLabel.label("Variant");
    private Label annotationLabel = DynamicLabel.label("Annotation");
    private Label symbolLabel = DynamicLabel.label("Symbol");
    private Label canonicalLabel = DynamicLabel.label("Canonical");
    private Label featureLabel = DynamicLabel.label("Feature");
    private File neo4jDBPath;
    private GraphDatabaseService graphDb;
    private VCFFileReader vcfFileReader;
    private ArrayList<VariantContext> variants = new ArrayList<>();
    private HashMap<String, Long> patientNodeIds = new HashMap<>(); //PatientID:NodeID
    private HashMap<String, Long> sampleNodeIds = new HashMap<>(); //sampleID:NodeID
    private HashMap<String, Long> variantNodeIds = new HashMap<>(); //Variant:NodeID
    private HashMap<String, Long> annotationNodeIds = new HashMap<>(); //AnnotationID:NodeID
    private HashMap<String, Long> geneNodeIds = new HashMap<>(); //GeneID:NodeID
    private HashMap<String, Long> featureNodeIds = new HashMap<>(); //Feature:NodeID
    private HashMap<String, HashSet<VEPAnnotation>> functionalAnnotations = new HashMap<>(); //variant:@annotation

    public VariantDatabase(VCFFileReader vcfFileReader, File neo4jDBPath){
        this.vcfFileReader = vcfFileReader;
        this.neo4jDBPath = neo4jDBPath;
    }

    private enum relTypes implements RelationshipType
    {
        HAS_SAMPLE,
        HAS_HET_VARIANT,
        HAS_HOM_VARIANT,
        IN_SYMBOL,
        IN_FEATURE,
        HAS_TRANSCRIPT_ABLATION_CONSEQUENCE,
        HAS_SPLICE_ACCEPTOR_VARIANT_CONSEQUENCE,
        HAS_SPLICE_DONOR_VARIANT_CONSEQUENCE,
        HAS_STOP_GAINED_CONSEQUENCE,
        HAS_FRAMESHIFT_VARIANT_CONSEQUENCE,
        HAS_STOP_LOST_CONSEQUENCE,
        HAS_START_LOST_CONSEQUENCE,
        HAS_TRANSCRIPT_AMPLIFICATION_CONSEQUENCE,
        HAS_INFRAME_INSERTION_CONSEQUENCE,
        HAS_INFRAME_DELETION_CONSEQUENCE,
        HAS_MISSENSE_VARIANT_CONSEQUENCE,
        HAS_PROTEIN_ALTERING_VARIANT_CONSEQUENCE,
        HAS_SPLICE_REGION_VARIANT_CONSEQUENCE,
        HAS_INCOMPLETE_TERMINAL_CODON_VARIANT_CONSEQUENCE,
        HAS_STOP_RETAINED_VARIANT_CONSEQUENCE,
        HAS_SYNONYMOUS_VARIANT_CONSEQUENCE,
        HAS_CODING_SEQUENCE_VARIANT_CONSEQUENCE,
        HAS_MATURE_MIRNA_VARIANT_CONSEQUENCE,
        HAS_5_PRIME_UTR_VARIANT_CONSEQUENCE,
        HAS_3_PRIME_UTR_VARIANT_CONSEQUENCE,
        HAS_NON_CODING_TRANSCRIPT_EXON_VARIANT_CONSEQUENCE,
        HAS_INTRON_VARIANT_CONSEQUENCE,
        HAS_NMD_TRANSCRIPT_VARIANT_CONSEQUENCE,
        HAS_NON_CODING_TRANSCRIPT_VARIANT_CONSEQUENCE,
        HAS_UPSTREAM_GENE_VARIANT_CONSEQUENCE,
        HAS_DOWNSTREAM_GENE_VARIANT_CONSEQUENCE,
        HAS_TFBS_ABLATION_CONSEQUENCE,
        HAS_TFBS_AMPLIFICATION_CONSEQUENCE,
        HAS_TF_BINDING_SITE_VARIANT_CONSEQUENCE,
        HAS_REGULATORY_REGION_ABLATION_CONSEQUENCE,
        HAS_REGULATORY_REGION_AMPLIFICATION_CONSEQUENCE,
        HAS_FEATURE_ELONGATION_CONSEQUENCE,
        HAS_REGULATORY_REGION_VARIANT_CONSEQUENCE,
        HAS_FEATURE_TRUNCATION_CONSEQUENCE,
        HAS_INTERGENIC_VARIANT_CONSEQUENCE
    }

    public void createDatabase(){

        //TODO check this is OK for opening existing DBs

        log.log(Level.INFO, "Starting database ...");

        //create DB
        graphDb = new GraphDatabaseFactory()
                .newEmbeddedDatabaseBuilder(neo4jDBPath)
                .newGraphDatabase();

        Neo4j.registerShutdownHook(graphDb);
    }

    public void createIndexes() {
        Neo4j.createConstraint(graphDb, patientLabel, "PatientID");
        Neo4j.createConstraint(graphDb, sampleLabel, "SampleID");
        Neo4j.createConstraint(graphDb, variantLabel, "Variant");
        Neo4j.createConstraint(graphDb, symbolLabel, "Symbol");
        Neo4j.createConstraint(graphDb, featureLabel, "Feature");
    }

    public void loadVCF() throws InvalidPropertiesFormatException {

        int n = 0;
        String variantLookup = "";

        log.log(Level.INFO, "Loading VCF into memory ...");

        Iterator<VariantContext> variantContextIterator = vcfFileReader.iterator();

        //read VCF file
        while (variantContextIterator.hasNext()) {

            n++;

            if (n > 250){
                break;
            }

            VariantContext variant = variantContextIterator.next();

            //skip hom-ref & filtered
            if (variant.isFiltered() || !variant.isVariant()){
                continue;
            }

            System.out.println(variant.getReference() + "\t" + variant.getAlternateAlleles().toString());

            //bank variants
            variants.add(variant);

            //loop over genotypes
            Iterator<Genotype> genotypeIterator = variant.getGenotypes().iterator();

            while (genotypeIterator.hasNext()) {
                Genotype genotype = genotypeIterator.next();

                boolean hasAnnotation = false;

                //skip wildtype or no calls
                if (genotype.isNoCall() || genotype.isHomRef()){
                    continue;
                }
                if (genotype.getPloidy() != 2) {
                    throw new InvalidPropertiesFormatException("Allele " + genotype.getAlleles().toString() + " is not diploid");
                }

                if(genotype.isHom()){
                    variantLookup = variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + genotype.getAllele(1).getBaseString();
                } else{
                    variantLookup = variant.getContig() + ":" + variant.getStart() + genotype.getAllele(0).getBaseString() + ">" + genotype.getAllele(1).getBaseString();
                }

                //initalise hash
                if (!functionalAnnotations.containsKey(variantLookup)){
                    functionalAnnotations.put(variantLookup, new HashSet<VEPAnnotation>());
                }

                //split annotations and make unique; group by alternative allele
                try {

                    //one annotation
                    VEPAnnotation vepAnnotation = new VEPAnnotation((String) variant.getAttribute("CSQ"));
                    vepAnnotation.parseAnnotation();

                    //transcript, ensg, hgnc
                    if (vepAnnotation.getFeature() == null || vepAnnotation.getGene() == null || vepAnnotation.getSymbol() == null) {
                        continue;
                    }

                    if(genotype.getAllele(0).getBaseString().length() == 1 && genotype.getAllele(1).getBaseString().length() == 1){
                        if (vepAnnotation.getAllele().equals(genotype.getAllele(1).getBaseString())){
                            functionalAnnotations.get(variantLookup).add(vepAnnotation);
                        }
                    } else if (genotype.getAllele(0).getBaseString().length() > 1 && genotype.getAllele(1).getBaseString().length() == 1){ //deletion
                        if (vepAnnotation.getAllele().equals("-")){
                            functionalAnnotations.get(variantLookup).add(vepAnnotation);
                        }
                    } else if (genotype.getAllele(0).getBaseString().length() == 1 && genotype.getAllele(1).getBaseString().length() > 1){ //insertion
                        if (vepAnnotation.getAllele().equals(genotype.getAllele(1).getBaseString().substring(1))){
                            functionalAnnotations.get(variantLookup).add(vepAnnotation);
                        }
                    } else {
                        throw new InvalidPropertiesFormatException("Not sure what is variant is: " + variantLookup);
                    }

                    hasAnnotation = true;

                } catch (ClassCastException e) {

                    //multiple annotations
                    for (String annotation : (ArrayList<String>) variant.getAttribute("CSQ")) {

                        //multiple annotations
                        VEPAnnotation vepAnnotation = new VEPAnnotation(annotation);
                        vepAnnotation.parseAnnotation();

                        //transcript, ensg, hgnc
                        if (vepAnnotation.getFeature() == null || vepAnnotation.getGene() == null || vepAnnotation.getSymbol() == null) {
                            continue;
                        }

                        if(genotype.getAllele(0).getBaseString().length() == 1 && genotype.getAllele(1).getBaseString().length() == 1){
                            if (vepAnnotation.getAllele().equals(genotype.getAllele(1).getBaseString())){
                                functionalAnnotations.get(variantLookup).add(vepAnnotation);
                            }
                        } else if (genotype.getAllele(0).getBaseString().length() > 1 && genotype.getAllele(1).getBaseString().length() == 1){ //deletion
                            if (vepAnnotation.getAllele().equals("-")){
                                functionalAnnotations.get(variantLookup).add(vepAnnotation);
                            }
                        } else if (genotype.getAllele(0).getBaseString().length() == 1 && genotype.getAllele(1).getBaseString().length() > 1){ //insertion
                            if (vepAnnotation.getAllele().equals(genotype.getAllele(1).getBaseString().substring(1))){
                                functionalAnnotations.get(variantLookup).add(vepAnnotation);
                            }
                        } else {
                            throw new InvalidPropertiesFormatException("Not sure what is variant is: " + variantLookup);
                        }

                        hasAnnotation = true;

                    }

                }

                //register missing alleles
                if (!hasAnnotation) {
                    log.log(Level.WARNING, "Missing annotation for: " + variantLookup); //TODO mop up and discprepent vars
                }

            }
        }
    }

    public void addPatientNodes(){
        log.log(Level.INFO, "Adding patients nodes ...");
        //TODO add patient nodes and link to samples
    }

    public void addSampleNodes(){

        log.log(Level.INFO, "Adding sample nodes ...");

        for (String sampleID : vcfFileReader.getFileHeader().getSampleNamesInOrder()){

            ArrayList<Node> nodes = Neo4j.getNodes(graphDb, sampleLabel, "SampleID", sampleID);

            if (nodes.size() == 0){

                //add sampleID
                try (Transaction tx = graphDb.beginTx()) {

                    Node node = graphDb.createNode();
                    node.addLabel(sampleLabel);
                    node.setProperty("SampleID", sampleID);

                    //store node Id for later
                    sampleNodeIds.put(sampleID, node.getId());

                    tx.success();
                }

            } else {
                sampleNodeIds.put(sampleID, nodes.get(0).getId());
            }

        }

    }

    public void addVariantsNodes() throws InvalidPropertiesFormatException {

        log.log(Level.INFO, "Adding variants ...");

        //read VCF records
        for (VariantContext variant : variants) {
            Iterator<Genotype> genotypeIterator = variant.getGenotypes().iterator();

            while (genotypeIterator.hasNext()) {
                Genotype genotype = genotypeIterator.next();

                //skip wildtype or no calls
                if (genotype.isNoCall() || genotype.isHomRef()) continue;

                String variantLookup = variant.getContig() + ":" + variant.getStart() + genotype.getAllele(0).getBaseString() + ">" + genotype.getAllele(1).getBaseString();
                ArrayList<Node> nodes = Neo4j.getNodes(graphDb, variantLabel, "Variant", variantLookup);

                if (nodes.size() == 0) {

                    //add variant
                    try (Transaction tx = graphDb.beginTx()) {

                        Node node = graphDb.createNode();
                        node.addLabel(variantLabel);

                        node.setProperty("Variant", variantLookup);

                        //store node Id for later
                        variantNodeIds.put(variantLookup, node.getId());

                        tx.success();
                    }

                } else {
                    //store node Id for later
                    variantNodeIds.put(variantLookup, nodes.get(0).getId());
                }
            }
        }

    }

    public void addGenotypeRelationships() throws InvalidPropertiesFormatException {

        log.log(Level.INFO, "Adding genotypes ...");

        String runID = "150716_D00501_0047_BHB092ADXX"; //TODO get library from VCF

        //read VCF records
        for (VariantContext variant : variants) {
            Iterator<Genotype> genotypeIterator = variant.getGenotypes().iterator();

            while (genotypeIterator.hasNext()) {
                Genotype genotype = genotypeIterator.next();

                //skip wildtype or no calls
                if (genotype.isNoCall() || genotype.isHomRef()) continue;

                String variantLookup = variant.getContig() + ":" + variant.getStart() + genotype.getAllele(0).getBaseString() + ">" + genotype.getAllele(1).getBaseString();

                //Add het variant
                if (genotype.isHet()){
                    try (Transaction tx = graphDb.beginTx()) {

                        Node node1 = graphDb.getNodeById(sampleNodeIds.get(genotype.getSampleName()));
                        Node node2 = graphDb.getNodeById(variantNodeIds.get(variantLookup));

                        Relationship relationship = node1.createRelationshipTo(node2, relTypes.HAS_HET_VARIANT);
                        relationship.setProperty("GQ", genotype.getGQ());
                        relationship.setProperty("RunID", runID);

                        tx.success();
                    }
                }

                //Add hom variant
                if (genotype.isHomVar()){
                    try (Transaction tx = graphDb.beginTx()) {

                        Node node1 = graphDb.getNodeById(sampleNodeIds.get(genotype.getSampleName()));
                        Node node2 = graphDb.getNodeById(variantNodeIds.get(variantLookup));

                        Relationship relationship = node1.createRelationshipTo(node2, relTypes.HAS_HOM_VARIANT);
                        relationship.setProperty("GQ", genotype.getGQ());
                        relationship.setProperty("RunID", runID);

                        tx.success();
                    }
                }

                //throw mixed
                if (genotype.isMixed()){
                    throw new InvalidPropertiesFormatException("Allele " + genotype.getAlleles().toString() + " is mixed");
                }

            }
        }

    }

    public void addSymbolNodes() throws InvalidPropertiesFormatException {

        log.log(Level.INFO, "Adding symbols ...");
        ArrayList<Node> nodes;

        //loop over variants
        for (VariantContext variant : variants) {

            //loop over genotypes
            Iterator<Genotype> genotypeIterator = variant.getGenotypes().iterator();

            while (genotypeIterator.hasNext()) {
                Genotype genotype = genotypeIterator.next();

                //skip wildtype or no calls
                if (genotype.isNoCall() || genotype.isHomRef()) continue;

                String variantLookup = variant.getContig() + ":" + variant.getStart() + genotype.getAllele(0).getBaseString() + ">" + genotype.getAllele(1).getBaseString();

                //loop over all annotations by for this allele
                for (VEPAnnotation annotation : functionalAnnotations.get(variantLookup)) {

                    //add symbol node
                    nodes = Neo4j.getNodes(graphDb, symbolLabel, "Symbol", annotation.getSymbol());
                    if (nodes.size() == 0){

                        try (Transaction tx = graphDb.beginTx()) {

                            Node node = graphDb.createNode();
                            node.addLabel(symbolLabel);

                            node.setProperty("GeneID", annotation.getGene());
                            node.setProperty("Symbol", annotation.getSymbol());
                            if (annotation.getSymbolSource() != null) node.setProperty("SymbolSource", annotation.getSymbolSource());

                            geneNodeIds.put(annotation.getSymbol(), node.getId());

                            tx.success();
                        }

                    } else {
                        geneNodeIds.put(annotation.getSymbol(), nodes.get(0).getId());
                    }

                }
            }
        }

    }

    public void addFeatureNodes() throws InvalidPropertiesFormatException {

        log.log(Level.INFO, "Adding features ...");
        ArrayList<Node> nodes;

        //loop over variants
        for (VariantContext variant : variants) {

            //loop over genotypes
            Iterator<Genotype> genotypeIterator = variant.getGenotypes().iterator();

            while (genotypeIterator.hasNext()) {
                Genotype genotype = genotypeIterator.next();

                //skip wildtype or no calls
                if (genotype.isNoCall() || genotype.isHomRef()) continue;

                String variantLookup = variant.getContig() + ":" + variant.getStart() + genotype.getAllele(0).getBaseString() + ">" + genotype.getAllele(1).getBaseString();

                //add transcript node
                for (VEPAnnotation annotation : functionalAnnotations.get(variantLookup)) {

                    nodes = Neo4j.getNodes(graphDb, featureLabel, "Feature", annotation.getFeature());
                    if (nodes.size() == 0) {

                        try (Transaction tx = graphDb.beginTx()) {

                            Node node = graphDb.createNode();
                            node.addLabel(featureLabel);

                            node.setProperty("Feature", annotation.getFeature());
                            if (annotation.getFeatureType() != null)
                                node.setProperty("FeatureType", annotation.getFeatureType());
                            if (annotation.getBiotype() != null) node.setProperty("Biotype", annotation.getBiotype());
                            if (annotation.getCanonical() != null && annotation.getCanonical().equals("YES"))
                                node.addLabel(canonicalLabel);

                            featureNodeIds.put(annotation.getFeature(), node.getId());

                            tx.success();
                        }

                    } else {
                        featureNodeIds.put(annotation.getFeature(), nodes.get(0).getId());
                    }
                }

            }
        }

    }

    public void addFunctionalAnnotationNodes() throws InvalidPropertiesFormatException {

        log.log(Level.INFO, "Adding functional annotations ...");
        ArrayList<Node> nodes;

        //loop over variants
        for (VariantContext variant : variants){

            //loop over genotypes
            Iterator<Genotype> genotypeIterator = variant.getGenotypes().iterator();

            while (genotypeIterator.hasNext()) {
                Genotype genotype = genotypeIterator.next();

                //skip wildtype or no calls
                if (genotype.isNoCall() || genotype.isHomRef()) continue;

                String variantLookup = variant.getContig() + ":" + variant.getStart() + genotype.getAllele(0).getBaseString() + ">" + genotype.getAllele(1).getBaseString();

                //loop over all annotations by for this allele
                for (VEPAnnotation annotation : functionalAnnotations.get(variantLookup)){

                    //add annotation node
                    nodes = getAnnotationsForVariant(variantLookup, graphDb);
                    if (nodes.size() == 0){

                        try (Transaction tx = graphDb.beginTx()) {

                            Node node = graphDb.createNode();
                            node.addLabel(annotationLabel);

                            if(annotation.getExon() != null) {
                                String[] fields = annotation.getExon().split("/");
                                node.setProperty("Exon", fields[0]);
                                node.setProperty("TotalExons", fields[1]);
                            }
                            if(annotation.getHgvsCoding() != null) node.setProperty("HGVSc", annotation.getHgvsCoding());
                            if(annotation.getHgvsProtein() != null) node.setProperty("HGVSp", annotation.getHgvsProtein());
                            if(annotation.getIntron() != null) {
                                String[] fields = annotation.getIntron().split("/");
                                node.setProperty("Intron", fields[0]);
                                node.setProperty("TotalIntrons", fields[1]);
                            }
                            if(annotation.getStrand() != null) node.setProperty("Strand", annotation.getStrand());
                            if(annotation.getPolyphen() != null) node.setProperty("PolyPhen", annotation.getPolyphen());
                            if(annotation.getSift() != null) node.setProperty("SIFT", annotation.getSift());

                            annotationNodeIds.put(variantLookup + ":" + annotation.getFeature(), node.getId());

                            tx.success();
                        }

                    } else {
                        if (nodes.size() == 1){
                            annotationNodeIds.put(variantLookup + ":" + annotation.getFeature(), nodes.get(0).getId());
                        } else {
                            throw new InvalidPropertiesFormatException("Variant " + variantLookup + " has non-unique annotations.");
                        }
                    }

                } //done looping over annotations

            } //done looping over genotypes

        }

    }

    public void addConsequenceRelationships() throws InvalidPropertiesFormatException {

        log.log(Level.INFO, "Linking Variants To Functional Annotations ...");

        //loop over variants
        for (VariantContext variant : variants) {

            //loop over genotypes
            Iterator<Genotype> genotypeIterator = variant.getGenotypes().iterator();

            while (genotypeIterator.hasNext()) {
                Genotype genotype = genotypeIterator.next();

                //skip wildtype or no calls
                if (genotype.isNoCall() || genotype.isHomRef()) continue;

                String variantLookup = variant.getContig() + ":" + variant.getStart() + genotype.getAllele(0).getBaseString() + ">" + genotype.getAllele(1).getBaseString();

                //loop over all annotations by for this allele
                for (VEPAnnotation annotation : functionalAnnotations.get(variantLookup)){

                    if (annotation.getConsequences().contains("transcript_ablation")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_TRANSCRIPT_ABLATION_CONSEQUENCE);
                    if (annotation.getConsequences().contains("splice_acceptor_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_SPLICE_ACCEPTOR_VARIANT_CONSEQUENCE);
                    if (annotation.getConsequences().contains("splice_donor_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_SPLICE_DONOR_VARIANT_CONSEQUENCE);
                    if (annotation.getConsequences().contains("stop_gained")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_STOP_GAINED_CONSEQUENCE);
                    if (annotation.getConsequences().contains("frameshift_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_FRAMESHIFT_VARIANT_CONSEQUENCE);
                    if (annotation.getConsequences().contains("stop_lost")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_STOP_LOST_CONSEQUENCE);
                    if (annotation.getConsequences().contains("start_lost")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_START_LOST_CONSEQUENCE);
                    if (annotation.getConsequences().contains("transcript_amplification")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_TRANSCRIPT_AMPLIFICATION_CONSEQUENCE);
                    if (annotation.getConsequences().contains("inframe_insertion")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_INFRAME_INSERTION_CONSEQUENCE);
                    if (annotation.getConsequences().contains("inframe_deletion")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_INFRAME_DELETION_CONSEQUENCE);
                    if (annotation.getConsequences().contains("missense_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_MISSENSE_VARIANT_CONSEQUENCE);
                    if (annotation.getConsequences().contains("protein_altering_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_PROTEIN_ALTERING_VARIANT_CONSEQUENCE);
                    if (annotation.getConsequences().contains("splice_region_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_SPLICE_REGION_VARIANT_CONSEQUENCE);
                    if (annotation.getConsequences().contains("incomplete_terminal_codon_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_INCOMPLETE_TERMINAL_CODON_VARIANT_CONSEQUENCE);
                    if (annotation.getConsequences().contains("stop_retained_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_STOP_RETAINED_VARIANT_CONSEQUENCE);
                    if (annotation.getConsequences().contains("synonymous_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_SYNONYMOUS_VARIANT_CONSEQUENCE);
                    if (annotation.getConsequences().contains("coding_sequence_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_CODING_SEQUENCE_VARIANT_CONSEQUENCE);
                    if (annotation.getConsequences().contains("mature_miRNA_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_MATURE_MIRNA_VARIANT_CONSEQUENCE);
                    if (annotation.getConsequences().contains("5_prime_UTR_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_5_PRIME_UTR_VARIANT_CONSEQUENCE);
                    if (annotation.getConsequences().contains("3_prime_UTR_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_3_PRIME_UTR_VARIANT_CONSEQUENCE);
                    if (annotation.getConsequences().contains("non_coding_transcript_exon_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_NON_CODING_TRANSCRIPT_EXON_VARIANT_CONSEQUENCE);
                    if (annotation.getConsequences().contains("intron_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_INTRON_VARIANT_CONSEQUENCE);
                    if (annotation.getConsequences().contains("NMD_transcript_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_NMD_TRANSCRIPT_VARIANT_CONSEQUENCE);
                    if (annotation.getConsequences().contains("non_coding_transcript_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_NON_CODING_TRANSCRIPT_VARIANT_CONSEQUENCE);
                    if (annotation.getConsequences().contains("upstream_gene_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_UPSTREAM_GENE_VARIANT_CONSEQUENCE);
                    if (annotation.getConsequences().contains("downstream_gene_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_DOWNSTREAM_GENE_VARIANT_CONSEQUENCE);
                    if (annotation.getConsequences().contains("TFBS_ablation")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_TFBS_ABLATION_CONSEQUENCE);
                    if (annotation.getConsequences().contains("TFBS_amplification")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_TFBS_AMPLIFICATION_CONSEQUENCE);
                    if (annotation.getConsequences().contains("TF_binding_site_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_TF_BINDING_SITE_VARIANT_CONSEQUENCE);
                    if (annotation.getConsequences().contains("regulatory_region_ablation")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_REGULATORY_REGION_ABLATION_CONSEQUENCE);
                    if (annotation.getConsequences().contains("regulatory_region_amplification")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_REGULATORY_REGION_AMPLIFICATION_CONSEQUENCE);
                    if (annotation.getConsequences().contains("feature_elongation")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_FEATURE_ELONGATION_CONSEQUENCE);
                    if (annotation.getConsequences().contains("regulatory_region_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_REGULATORY_REGION_VARIANT_CONSEQUENCE);
                    if (annotation.getConsequences().contains("feature_truncation")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_FEATURE_TRUNCATION_CONSEQUENCE);
                    if (annotation.getConsequences().contains("intergenic_variant")) Neo4j.createRelationship(graphDb, variantNodeIds.get(variantLookup), annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()), relTypes.HAS_INTERGENIC_VARIANT_CONSEQUENCE);

                }

            }
        }

    }

    public void addInSymbolRelationships() throws InvalidPropertiesFormatException {

        log.log(Level.INFO, "Linking features to symbol ...");

        //loop over variants
        for (VariantContext variant : variants) {

            //loop over genotypes
            Iterator<Genotype> genotypeIterator = variant.getGenotypes().iterator();

            while (genotypeIterator.hasNext()) {
                Genotype genotype = genotypeIterator.next();

                //skip wildtype or no calls
                if (genotype.isNoCall() || genotype.isHomRef()) continue;

                String variantLookup = variant.getContig() + ":" + variant.getStart() + genotype.getAllele(0).getBaseString() + ">" + genotype.getAllele(1).getBaseString();

                //loop over all annotations by for this variant
                for (VEPAnnotation annotation : functionalAnnotations.get(variantLookup)) {

                    //link feature to symbol
                    try (Transaction tx = graphDb.beginTx()) {

                        Node node1 = graphDb.getNodeById(featureNodeIds.get(annotation.getFeature()));
                        Node node2 = graphDb.getNodeById(geneNodeIds.get(annotation.getSymbol()));

                        Relationship relationship = node1.createRelationshipTo(node2, relTypes.IN_SYMBOL);

                        tx.success();
                    }

                }

            }
        }

    }

    public void addInFeatureRelationships() throws InvalidPropertiesFormatException {

        log.log(Level.INFO, "Linking annotations to features ...");

        ArrayList<Node> nodes;

        //loop over variants
        for (VariantContext variant : variants) {

            //loop over genotypes
            Iterator<Genotype> genotypeIterator = variant.getGenotypes().iterator();

            while (genotypeIterator.hasNext()) {
                Genotype genotype = genotypeIterator.next();

                //skip wildtype or no calls
                if (genotype.isNoCall() || genotype.isHomRef()) continue;

                String variantLookup = variant.getContig() + ":" + variant.getStart() + genotype.getAllele(0).getBaseString() + ">" + genotype.getAllele(1).getBaseString();

                //loop over all annotations by for this allele
                for (VEPAnnotation annotation : functionalAnnotations.get(variantLookup)) {

                    //add in feature relationship
                    nodes = getSymbolsForFeature(annotation.getFeature(), graphDb);
                    if (nodes.size() == 0){

                        //link annotation to feature
                        try (Transaction tx = graphDb.beginTx()) {

                            Node node1 = graphDb.getNodeById(annotationNodeIds.get(variantLookup + ":" + annotation.getFeature()));
                            Node node2 = graphDb.getNodeById(featureNodeIds.get(annotation.getFeature()));

                            Relationship relationship = node1.createRelationshipTo(node2, relTypes.IN_FEATURE);
                            tx.success();

                        }

                    } else if (nodes.size() > 1) {
                        throw new InvalidPropertiesFormatException("Variant " + variantLookup + " has non-unique annotations.");
                    }

                }
            }
        }

    }

    private static void convertToMinimalRepresentation(){
        //TODO
    }

    private static ArrayList<Node> getAnnotationsForVariant(String variantLookup, final GraphDatabaseService graphDb){
        ArrayList<Node> nodes = new ArrayList<>();

        for (Map<String, Object> result : Neo4j.runCypherQuery(graphDb,
                "MATCH (v:Variant {Variant:\"" + variantLookup + "\"})-[]-(a:Annotation) return distinct a;"
        )){
            for (Map.Entry<String, Object> iter : result.entrySet()){
                nodes.add((Node) iter.getValue());
            }
        }

        return nodes;
    }

    private static ArrayList<Node> getSymbolsForFeature(String feature, final GraphDatabaseService graphDb){
        ArrayList<Node> nodes = new ArrayList<>();

        for (Map<String, Object> result : Neo4j.runCypherQuery(graphDb,
                "MATCH (f:Feature {Feature:\"" + feature + "\"})-[]-(s:Symbol) return distinct s;"
        )){
            for (Map.Entry<String, Object> iter : result.entrySet()){
                nodes.add((Node) iter.getValue());
            }
        }

        return nodes;
    }

    public void shutdownDatabase(){
        log.log(Level.INFO, "Shutting down database ...");
        graphDb.shutdown();
    }

}