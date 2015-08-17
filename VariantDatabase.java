package nhs.genetics.cardiff;

import htsjdk.variant.variantcontext.Allele;
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
    private Label featureLabel = DynamicLabel.label("Feature");
    private File neo4jDBPath;
    private GraphDatabaseService graphDb;
    private VCFFileReader vcfFileReader;
    private ArrayList<VariantContext> variants = new ArrayList<>();
    private HashMap<String, Long> patientNodeIds = new HashMap<>(); //PatientID:NodeID
    private HashMap<String, Long> sampleNodeIds = new HashMap<>(); //sampleID:NodeID
    private HashMap<String, Long> variantNodeIds = new HashMap<>(); //VariantID:NodeID
    private HashMap<String, Long> annotationNodeIds = new HashMap<>(); //AnnotationID:NodeID
    private HashMap<String, Long> geneNodeIds = new HashMap<>(); //GeneID:NodeID
    private HashMap<String, Long> featureNodeIds = new HashMap<>(); //Feature:NodeID

    public VariantDatabase(VCFFileReader vcfFileReader, File neo4jDBPath){
        this.vcfFileReader = vcfFileReader;
        this.neo4jDBPath = neo4jDBPath;
    }

    private enum relTypes implements RelationshipType
    {
        HAS_SAMPLE,
        HAS_HET_VARIANT,
        HAS_HETNONREF_VARIANT,
        HAS_HOMVAR_VARIANT,
        HAS_MIXED_VARIANT,
        IN_SYMBOL,
        HAS_FEATURE,
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
        Neo4j.createConstraint(graphDb, variantLabel, "VariantID");
        Neo4j.createConstraint(graphDb, annotationLabel, "AnnotationID");
        Neo4j.createConstraint(graphDb, symbolLabel, "Symbol");
        Neo4j.createConstraint(graphDb, featureLabel, "Feature");
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

    public void loadVariantsIntoMemory(){

        int n = 0;

        log.log(Level.INFO, "Loading variants into memory ...");

        Iterator<VariantContext> variantContextIterator = vcfFileReader.iterator();

        //read VCF file
        while (variantContextIterator.hasNext()) {
            VariantContext variant = variantContextIterator.next();
            n++;

            if (n > 1000){
                break;
            }

            //skip hom-ref & filtered
            if (!variant.isFiltered() && variant.isVariant()){
                variants.add(variant);
            }

        }


    }

    //TODO partition by RunID
    public void addVariantsAndLinkSamples(){

        log.log(Level.INFO, "Adding variants and linking samples ...");

        //read VCF file
        for (VariantContext variant : variants){
            Iterator<Genotype> genotypeIterator = variant.getGenotypes().iterator();

            while (genotypeIterator.hasNext()){
                Genotype genotype = genotypeIterator.next();

                //skip wildtype or no calls
                if (genotype.isNoCall() || genotype.isHomRef()) continue;

                String variantLookup = variant.getContig() + ":" + variant.getStart() + genotype.getAllele(0).getBaseString() + ">" + genotype.getAllele(1).getBaseString();
                ArrayList<Node> nodes = Neo4j.getNodes(graphDb, variantLabel, "VariantID", variantLookup);

                if (nodes.size() == 0){

                    //add variant
                    try (Transaction tx = graphDb.beginTx()) {

                        Node node = graphDb.createNode();
                        node.addLabel(variantLabel);

                        node.setProperty("VariantID", variantLookup);
                        node.setProperty("Contig", variant.getContig());
                        node.setProperty("Position", variant.getStart());
                        node.setProperty("Reference", genotype.getAllele(0).getBaseString());
                        node.setProperty("Alternative", genotype.getAllele(1).getBaseString());

                        //store node Id for later
                        variantNodeIds.put(variantLookup, node.getId());

                        tx.success();
                    }

                } else {
                    //store node Id for later
                    variantNodeIds.put(variantLookup, nodes.get(0).getId());
                }

                //Add het variant
                if (genotype.isHet()){
                    try (Transaction tx = graphDb.beginTx()) {

                        Node node1 = graphDb.getNodeById(sampleNodeIds.get(genotype.getSampleName()));
                        Node node2 = graphDb.getNodeById(variantNodeIds.get(variantLookup));

                        Relationship relationship = node1.createRelationshipTo(node2, relTypes.HAS_HET_VARIANT);
                        relationship.setProperty("GQ", genotype.getGQ());

                        tx.success();
                    }
                }

                //Add hom variant
                if (genotype.isHet()){
                    try (Transaction tx = graphDb.beginTx()) {

                        Node node1 = graphDb.getNodeById(sampleNodeIds.get(genotype.getSampleName()));
                        Node node2 = graphDb.getNodeById(variantNodeIds.get(variantLookup));

                        Relationship relationship = node1.createRelationshipTo(node2, relTypes.HAS_HET_VARIANT);
                        relationship.setProperty("GQ", genotype.getGQ());

                        tx.success();
                    }
                }

                //Add mixed variant
                if (genotype.isHet()){
                    System.err.println("mixed var");
                    System.exit(1);
                }

            }

        }

    }

    //TODO use Ensembl v79 annotations, enable all VEP fields, add dbSNSFP
    //todo check consistency between link samples
    public void addFunctionalAnnotations(){

        log.log(Level.INFO, "Adding functional annotations ...");

        ArrayList<Node> nodes;
        HashMap<String, HashSet<VEPAnnotation>> splitAnnotations = new HashMap<>();

        //read VCF file
        for (VariantContext variant : variants){

            //split annotations and make unique
            try {

                //one annotation
                VEPAnnotation splitAnnotation = new VEPAnnotation((String)variant.getAttribute("CSQ"));
                splitAnnotation.parseAnnotation();

                if (!splitAnnotations.containsKey(splitAnnotation.getAllele())) splitAnnotations.put(splitAnnotation.getAllele(), new HashSet<VEPAnnotation>());
                splitAnnotations.get(splitAnnotation.getAllele()).add(splitAnnotation);

            } catch (ClassCastException e){

                //multiple annotations
                for (String annotation : (ArrayList<String>) variant.getAttribute("CSQ")){

                    //multiple annotations
                    VEPAnnotation splitAnnotation = new VEPAnnotation(annotation);
                    splitAnnotation.parseAnnotation();

                    if (!splitAnnotations.containsKey(splitAnnotation.getAllele())) splitAnnotations.put(splitAnnotation.getAllele(), new HashSet<VEPAnnotation>());
                    splitAnnotations.get(splitAnnotation.getAllele()).add(splitAnnotation);

                }

            }

            //loop over alleles
            for (Allele allele : variant.getAlleles()){

                //skip wildtype alleles
                if (allele.isReference()){
                    continue;
                }

                String variantLookup = variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString();

                //skip missing alleles
                if (!splitAnnotations.containsKey(allele.getBaseString())) {
                    log.log(Level.WARNING, "Missing annotation for: " + variantLookup);
                    continue;
                }

                //loop over annotations by for this allele
                for (VEPAnnotation annotation : splitAnnotations.get(allele.getBaseString())){

                    //transcript, ensg, hgnc
                    if (annotation.getFeature() == null || annotation.getGene() == null || annotation.getSymbol() == null){
                        continue;
                    }

                    //add gene node
                    nodes = Neo4j.getNodes(graphDb, symbolLabel, "Symbol", annotation.getSymbol());
                    if (nodes.size() == 0){

                        try (Transaction tx = graphDb.beginTx()) {

                            Node node = graphDb.createNode();
                            node.addLabel(symbolLabel);

                            node.setProperty("GeneID", annotation.getGene());
                            node.setProperty("Symbol", annotation.getSymbol());

                            geneNodeIds.put(annotation.getGene(), node.getId());

                            tx.success();
                        }

                    } else {
                        geneNodeIds.put(annotation.getSymbol(), nodes.get(0).getId());
                    }

                    //add transcript node
                    nodes = Neo4j.getNodes(graphDb, featureLabel, "Feature", annotation.getFeature());
                    if (nodes.size() == 0){

                        try (Transaction tx = graphDb.beginTx()) {

                            Node node = graphDb.createNode();
                            node.addLabel(featureLabel);

                            node.setProperty("Feature", annotation.getFeature());
                            node.setProperty("FeatureType", annotation.getFeatureType());

                            featureNodeIds.put(annotation.getFeature(), node.getId());

                            tx.success();
                        }

                    } else {
                        featureNodeIds.put(annotation.getFeature(), nodes.get(0).getId());
                    }

                    //add annotation node
                    //TODO check for multiple annotations for the same transcript not AnnotationID?
                    nodes = Neo4j.getNodes(graphDb, annotationLabel, "AnnotationID", variantLookup + ":" + annotation.getFeature());
                    if (nodes.size() == 0){

                        try (Transaction tx = graphDb.beginTx()) {

                            Node node = graphDb.createNode();
                            node.addLabel(annotationLabel);
                            node.setProperty("AnnotationID", variantLookup + ":" + annotation.getFeature());

                            if(annotation.getExon() != null) {
                                String[] fields = annotation.getExon().split("/");
                                node.setProperty("Exon", fields[0]);
                                node.setProperty("TotalExons", fields[1]);
                            }
                            if(annotation.getFeature() != null) node.setProperty("Feature", annotation.getFeature());
                            if(annotation.getFeatureType() != null) node.setProperty("FeatureType", annotation.getFeatureType());
                            if(annotation.getHgvsCoding() != null) node.setProperty("HGVSc", annotation.getHgvsCoding());
                            if(annotation.getHgvsProtein() != null) node.setProperty("HGVSp", annotation.getHgvsProtein());
                            if(annotation.getIntron() != null) {
                                String[] fields = annotation.getIntron().split("/");
                                node.setProperty("Intron", fields[0]);
                                node.setProperty("TotalIntrons", fields[1]);
                            }
                            if(annotation.getStrand() != null) node.setProperty("Strand", annotation.getStrand());

                            annotationNodeIds.put(variantLookup + ":" + annotation.getFeature(), node.getId());

                            tx.success();
                        }

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

                    } else {
                        annotationNodeIds.put(variantLookup + ":" + annotation.getFeature(),
                                Neo4j.getNodes(graphDb, annotationLabel, "AnnotationID", variantLookup + ":" + annotation.getFeature()).get(0).getId());
                    }

                } //done looping over annotations

            } //done looping over alleles

            splitAnnotations.clear();

        } //done reading VCF

    }

    private static void convertToMinimalRepresentation(){

    }

    public void shutdownDatabase(){
        log.log(Level.INFO, "Shutting down database ...");
        graphDb.shutdown();
    }

}
