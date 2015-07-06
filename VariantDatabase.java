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

    private Label sampleLabel = DynamicLabel.label("Sample");
    private Label variantLabel = DynamicLabel.label("Variant");
    private Label annotationLabel = DynamicLabel.label("Annotation");
    private Label geneLabel = DynamicLabel.label("Gene");
    private Label featureLabel = DynamicLabel.label("Feature");
    private Label patientLabel = DynamicLabel.label("Patient");

    private File neo4jDBPath;
    private GraphDatabaseService graphDb;
    private VCFFileReader vcfFileReader;
    private HashMap<String, Long> sampleNodeIds = new HashMap<>(); //sampleID:NodeID
    private HashMap<String, Long> variantNodeIds = new HashMap<>(); //VariantID:NodeID
    private HashMap<String, Long> annotationNodeIds = new HashMap<>(); //AnnotationID:NodeID
    private HashMap<String, Long> geneNodeIds = new HashMap<>(); //GeneID:NodeID
    private HashMap<String, Long> patientNodeIds = new HashMap<>(); //PatientID:NodeID
    private HashMap<String, Long> featureNodeIds = new HashMap<>(); //Feature:NodeID

    public VariantDatabase(VCFFileReader vcfFileReader, File neo4jDBPath){
        this.vcfFileReader = vcfFileReader;
        this.neo4jDBPath = neo4jDBPath;
    }

    private enum relTypes implements RelationshipType
    {
        HAS_SAMPLE,
        HAS_VARIANT,
        IN_GENE,
        HAS_FEATURE,
        HAS_ANNOTATION
    }

    public void createDatabase(){

        log.log(Level.INFO, "Starting database ...");

        //create DB
        graphDb = new GraphDatabaseFactory()
                .newEmbeddedDatabaseBuilder(neo4jDBPath)
                .newGraphDatabase();

        Neo4j.registerShutdownHook(graphDb);
    }

    public void createIndexes(){
        Neo4j.createConstraint(graphDb, sampleLabel, "SampleID");
        Neo4j.createConstraint(graphDb, variantLabel, "VariantID");
        Neo4j.createConstraint(graphDb, annotationLabel, "AnnotationID");
        Neo4j.createConstraint(graphDb, geneLabel, "Symbol");
        Neo4j.createConstraint(graphDb, featureLabel, "Feature");
        Neo4j.createConstraint(graphDb, patientLabel, "PatientID");
    }

    public void addPatientNodes(){
        log.log(Level.INFO, "Adding patients nodes ...");
        //TODO
    }

    public void addSampleNodes(){

        log.log(Level.INFO, "Adding sample nodes ...");

        Node node;

        for (String sampleID : vcfFileReader.getFileHeader().getSampleNamesInOrder()){

            if (!Neo4j.hasNode(graphDb, sampleLabel, "SampleID", sampleID)){

                //add sampleID
                try (Transaction tx = graphDb.beginTx()) {

                    node = graphDb.createNode();
                    node.addLabel(sampleLabel);
                    node.setProperty("SampleID", sampleID);

                    //store node Id for later
                    sampleNodeIds.put(sampleID, node.getId());

                    tx.success();
                }

            } else {
                log.log(Level.WARNING, sampleID + " already exists in database.");
                //TODO: delete sample and relationships?
            }

        }

    }

    public void addVariants(){

        log.log(Level.INFO, "Adding variants ...");

        Node node;
        Iterator<VariantContext> variantContextIterator = vcfFileReader.iterator();

        //read VCF file
        while (variantContextIterator.hasNext()) {
            VariantContext variant = variantContextIterator.next();

            //skip hom-ref & filtered
            if (variant.isFiltered() || !variant.isVariant()) continue;

            //loop over alleles
            for (Allele allele : variant.getAlleles()){

                //skip wildtype alleles
                if (allele.isReference()) continue;

                if (!Neo4j.hasNode(graphDb, variantLabel, "VariantID",  variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString())){

                    //add variant
                    try (Transaction tx = graphDb.beginTx()) {

                        node = graphDb.createNode();
                        node.addLabel(variantLabel);

                        node.setProperty("VariantID", variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString());
                        node.setProperty("Contig", variant.getContig());
                        node.setProperty("Position", variant.getStart());
                        node.setProperty("Reference", variant.getReference().getBaseString());
                        node.setProperty("Alternative", allele.getBaseString());

                        //store node Id for later
                        variantNodeIds.put(variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString(), node.getId());

                        tx.success();
                    }

                } else {
                    //store node Id for later
                    variantNodeIds.put(variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString(),
                            Neo4j.getUniqueMergedNode(graphDb, variantLabel, "VariantID", variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString()).getId());
                }

            }

        }

    }

    public void addFunctionalAnnotations(){

        log.log(Level.INFO, "Adding functional annotations ...");

        Node node;
        HashMap<String, HashSet<VEPAnnotation>> splitAnnotations = new HashMap<>();
        Iterator<VariantContext> variantContextIterator = vcfFileReader.iterator();

        //read VCF file
        while (variantContextIterator.hasNext()) {
            VariantContext variant = variantContextIterator.next();

            //skip hom-ref & filtered
            if (variant.isFiltered() || !variant.isVariant()) continue;

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
                if (allele.isReference()) continue;

                //skip missing alleles
                if (!splitAnnotations.containsKey(allele.getBaseString())) {
                    log.log(Level.WARNING, "Missing annotation for: " + variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString());
                    continue;
                }

                //loop over annotations by for this allele
                for (VEPAnnotation annotation : splitAnnotations.get(allele.getBaseString())){

                    //transcript, ensg, symbol
                    if (annotation.getFeature() == null || annotation.getGene() == null || annotation.getSymbol() == null) continue;

                    //add gene node
                    if (!Neo4j.hasNode(graphDb, geneLabel, "Symbol", annotation.getSymbol())){

                        try (Transaction tx = graphDb.beginTx()) {

                            node = graphDb.createNode();
                            node.addLabel(geneLabel);

                            node.setProperty("GeneID", annotation.getGene());
                            node.setProperty("Symbol", annotation.getSymbol());

                            geneNodeIds.put(annotation.getGene(), node.getId());

                            tx.success();
                        }

                    } else {
                        geneNodeIds.put(annotation.getSymbol(), Neo4j.getUniqueMergedNode(graphDb, geneLabel, "Symbol", annotation.getSymbol()).getId());
                    }

                    //add transcript node
                    if (!Neo4j.hasNode(graphDb, featureLabel, "Feature", annotation.getFeature())){

                        try (Transaction tx = graphDb.beginTx()) {

                            node = graphDb.createNode();
                            node.addLabel(featureLabel);

                            node.setProperty("Feature", annotation.getFeature());
                            node.setProperty("FeatureType", annotation.getFeatureType());

                            featureNodeIds.put(annotation.getFeature(), node.getId());

                            tx.success();
                        }

                    } else {
                        featureNodeIds.put(annotation.getFeature(), Neo4j.getUniqueMergedNode(graphDb, featureLabel, "Symbol", annotation.getSymbol()).getId());
                    }

                    //make annotation node
                    if (!Neo4j.hasNode(graphDb, annotationLabel, "AnnotationID", variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString() + ":" + annotation.getFeature())){

                        try (Transaction tx = graphDb.beginTx()) {

                            node = graphDb.createNode();
                            node.addLabel(annotationLabel);
                            node.setProperty("AnnotationID", variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString() + ":" + annotation.getFeature());

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

                            annotationNodeIds.put(variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString() + ":" + annotation.getFeature(), node.getId());

                            tx.success();
                        }

                        //link variant to annotation
                        try (Transaction tx = graphDb.beginTx()) {

                            Relationship relationship = graphDb.getNodeById(variantNodeIds.get(variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString()))
                                    .createRelationshipTo(node, relTypes.HAS_ANNOTATION);

                            //transcript_ablation
                            if (annotation.getConsequences().contains("transcript_ablation")){
                                relationship.setProperty("transcript_ablation", true);
                            } else {
                                relationship.setProperty("transcript_ablation", false);
                            }

                            //splice_donor_variant
                            if (annotation.getConsequences().contains("splice_donor_variant")){
                                relationship.setProperty("splice_donor_variant", true);
                            } else {
                                relationship.setProperty("splice_donor_variant", false);
                            }

                            //splice_acceptor_variant
                            if (annotation.getConsequences().contains("splice_acceptor_variant")){
                                relationship.setProperty("splice_acceptor_variant", true);
                            } else {
                                relationship.setProperty("splice_acceptor_variant", false);
                            }

                            //stop_gained
                            if (annotation.getConsequences().contains("stop_gained")){
                                relationship.setProperty("stop_gained", true);
                            } else {
                                relationship.setProperty("stop_gained", false);
                            }

                            //frameshift_variant
                            if (annotation.getConsequences().contains("frameshift_variant")){
                                relationship.setProperty("frameshift_variant", true);
                            } else {
                                relationship.setProperty("frameshift_variant", false);
                            }

                            //stop_lost
                            if (annotation.getConsequences().contains("stop_lost")){
                                relationship.setProperty("stop_lost", true);
                            } else {
                                relationship.setProperty("stop_lost", false);
                            }

                            //initiator_codon_variant
                            if (annotation.getConsequences().contains("initiator_codon_variant")){
                                relationship.setProperty("initiator_codon_variant", true);
                            } else {
                                relationship.setProperty("initiator_codon_variant", false);
                            }

                            //transcript_amplification
                            if (annotation.getConsequences().contains("transcript_amplification")){
                                relationship.setProperty("transcript_amplification", true);
                            } else {
                                relationship.setProperty("transcript_amplification", false);
                            }

                            //inframe_insertion
                            if (annotation.getConsequences().contains("inframe_insertion")){
                                relationship.setProperty("inframe_insertion", true);
                            } else {
                                relationship.setProperty("inframe_insertion", false);
                            }

                            //inframe_deletion
                            if (annotation.getConsequences().contains("inframe_deletion")){
                                relationship.setProperty("inframe_deletion", true);
                            } else {
                                relationship.setProperty("inframe_deletion", false);
                            }

                            //missense_variant
                            if (annotation.getConsequences().contains("missense_variant")){
                                relationship.setProperty("missense_variant", true);
                            } else {
                                relationship.setProperty("missense_variant", false);
                            }

                            //splice_region_variant
                            if (annotation.getConsequences().contains("splice_region_variant")){
                                relationship.setProperty("splice_region_variant", true);
                            } else {
                                relationship.setProperty("splice_region_variant", false);
                            }

                            //incomplete_terminal_codon_variant
                            if (annotation.getConsequences().contains("incomplete_terminal_codon_variant")){
                                relationship.setProperty("incomplete_terminal_codon_variant", true);
                            } else {
                                relationship.setProperty("incomplete_terminal_codon_variant", false);
                            }

                            //stop_retained_variant
                            if (annotation.getConsequences().contains("stop_retained_variant")){
                                relationship.setProperty("stop_retained_variant", true);
                            } else {
                                relationship.setProperty("stop_retained_variant", false);
                            }

                            //synonymous_variant
                            if (annotation.getConsequences().contains("synonymous_variant")){
                                relationship.setProperty("synonymous_variant", true);
                            } else {
                                relationship.setProperty("synonymous_variant", false);
                            }

                            //coding_sequence_variant
                            if (annotation.getConsequences().contains("coding_sequence_variant")){
                                relationship.setProperty("coding_sequence_variant", true);
                            } else {
                                relationship.setProperty("coding_sequence_variant", false);
                            }

                            //mature_miRNA_variant
                            if (annotation.getConsequences().contains("mature_miRNA_variant")){
                                relationship.setProperty("mature_miRNA_variant", true);
                            } else {
                                relationship.setProperty("mature_miRNA_variant", false);
                            }

                            //five_prime_UTR_variant
                            if (annotation.getConsequences().contains("five_prime_UTR_variant")){
                                relationship.setProperty("five_prime_UTR_variant", true);
                            } else {
                                relationship.setProperty("five_prime_UTR_variant", false);
                            }

                            //three_prime_UTR_variant
                            if (annotation.getConsequences().contains("three_prime_UTR_variant")){
                                relationship.setProperty("three_prime_UTR_variant", true);
                            } else {
                                relationship.setProperty("three_prime_UTR_variant", false);
                            }

                            //non_coding_transcript_exon_variant
                            if (annotation.getConsequences().contains("non_coding_transcript_exon_variant")){
                                relationship.setProperty("non_coding_transcript_exon_variant", true);
                            } else {
                                relationship.setProperty("non_coding_transcript_exon_variant", false);
                            }

                            //intron_variant
                            if (annotation.getConsequences().contains("intron_variant")){
                                relationship.setProperty("intron_variant", true);
                            } else {
                                relationship.setProperty("intron_variant", false);
                            }

                            //NMD_transcript_variant
                            if (annotation.getConsequences().contains("NMD_transcript_variant")){
                                relationship.setProperty("NMD_transcript_variant", true);
                            } else {
                                relationship.setProperty("NMD_transcript_variant", false);
                            }

                            //non_coding_transcript_variant
                            if (annotation.getConsequences().contains("non_coding_transcript_variant")){
                                relationship.setProperty("non_coding_transcript_variant", true);
                            } else {
                                relationship.setProperty("non_coding_transcript_variant", false);
                            }

                            //upstream_gene_variant
                            if (annotation.getConsequences().contains("upstream_gene_variant")){
                                relationship.setProperty("upstream_gene_variant", true);
                            } else {
                                relationship.setProperty("upstream_gene_variant", false);
                            }

                            //downstream_gene_variant
                            if (annotation.getConsequences().contains("downstream_gene_variant")){
                                relationship.setProperty("downstream_gene_variant", true);
                            } else {
                                relationship.setProperty("downstream_gene_variant", false);
                            }

                            //TFBS_ablation
                            if (annotation.getConsequences().contains("TFBS_ablation")){
                                relationship.setProperty("TFBS_ablation", true);
                            } else {
                                relationship.setProperty("TFBS_ablation", false);
                            }

                            //TFBS_amplification
                            if (annotation.getConsequences().contains("TFBS_amplification")){
                                relationship.setProperty("TFBS_amplification", true);
                            } else {
                                relationship.setProperty("TFBS_amplification", false);
                            }

                            //TF_binding_site_variant
                            if (annotation.getConsequences().contains("TF_binding_site_variant")){
                                relationship.setProperty("TF_binding_site_variant", true);
                            } else {
                                relationship.setProperty("TF_binding_site_variant", false);
                            }

                            //regulatory_region_ablation
                            if (annotation.getConsequences().contains("regulatory_region_ablation")){
                                relationship.setProperty("regulatory_region_ablation", true);
                            } else {
                                relationship.setProperty("regulatory_region_ablation", false);
                            }

                            //regulatory_region_amplification
                            if (annotation.getConsequences().contains("regulatory_region_amplification")){
                                relationship.setProperty("regulatory_region_amplification", true);
                            } else {
                                relationship.setProperty("regulatory_region_amplification", false);
                            }

                            //regulatory_region_variant
                            if (annotation.getConsequences().contains("regulatory_region_variant")){
                                relationship.setProperty("regulatory_region_variant", true);
                            } else {
                                relationship.setProperty("regulatory_region_variant", false);
                            }

                            //feature_elongation
                            if (annotation.getConsequences().contains("feature_elongation")){
                                relationship.setProperty("feature_elongation", true);
                            } else {
                                relationship.setProperty("feature_elongation", false);
                            }

                            //feature_truncation
                            if (annotation.getConsequences().contains("feature_truncation")){
                                relationship.setProperty("feature_truncation", true);
                            } else {
                                relationship.setProperty("feature_truncation", false);
                            }

                            //intergenic_variant
                            if (annotation.getConsequences().contains("intergenic_variant")){
                                relationship.setProperty("intergenic_variant", true);
                            } else {
                                relationship.setProperty("intergenic_variant", false);
                            }

                            tx.success();

                        } //done creating relationship

                    } else {
                        annotationNodeIds.put(variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString() + ":" + annotation.getFeature(),
                                Neo4j.getUniqueMergedNode(graphDb, annotationLabel, "AnnotationID", variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString() + ":" + annotation.getFeature()).getId());
                    }

                } //done looping over annotations

            } //done looping over alleles

            splitAnnotations.clear();

        } //done reading VCF

    }

    //TODO check why? duplicate relationships
    public void linkSamplesToVariants(){

        log.log(Level.INFO, "Linking samples to variants ...");

        Iterator<VariantContext> variantContextIterator = vcfFileReader.iterator();

        //read VCF file
        while (variantContextIterator.hasNext()) {
            VariantContext variant = variantContextIterator.next();

            //skip hom-ref & filtered
            if (variant.isFiltered() || !variant.isVariant()) continue;

            //loop over genotypes
            Iterator<Genotype> genotypeIterator = variant.getGenotypes().iterator();

            while (genotypeIterator.hasNext()){
                Genotype genotype = genotypeIterator.next();

                //skip wildtype or no calls
                if (genotype.isNoCall() || genotype.isHomRef()) continue;

                //loop over alleles
                for (Allele allele : genotype.getAlleles()){

                    //skip wildtype alleles
                    if (allele.isReference()) continue;

                    //create relationship
                    try (Transaction tx = graphDb.beginTx()) {

                        Relationship relationship = graphDb.getNodeById(sampleNodeIds.get(genotype.getSampleName()))
                                .createRelationshipTo(graphDb.getNodeById(variantNodeIds.get(variant.getContig() + ":" + variant.getStart() + variant.getReference().getBaseString() + ">" + allele.getBaseString())), relTypes.HAS_VARIANT);

                        if (genotype.isHet()) relationship.setProperty("Het", true); else relationship.setProperty("Het", false);
                        if (genotype.isHetNonRef()) relationship.setProperty("HetNonRef", true); else relationship.setProperty("HetNonRef", false);
                        if (genotype.isHomVar()) relationship.setProperty("HomVar", true); else relationship.setProperty("HomVar", false);
                        if (genotype.isMixed()) relationship.setProperty("Mixed", true); else relationship.setProperty("Mixed", false);
                        relationship.setProperty("GQ", genotype.getGQ());

                        tx.success();

                    }

                }
            }

        }

    }

    public void shutdownDatabase(){
        log.log(Level.INFO, "Shutting down database ...");
        graphDb.shutdown();
    }

}
