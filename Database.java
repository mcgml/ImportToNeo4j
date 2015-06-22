package nhs.genetics.cardiff;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.neo4j.graphdb.*;
import org.neo4j.graphdb.schema.IndexDefinition;
import org.neo4j.graphdb.schema.Schema;
import org.neo4j.graphdb.Direction;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.RelationshipType;
import org.neo4j.graphdb.Transaction;
import org.neo4j.graphdb.traversal.Evaluators;
import org.neo4j.graphdb.traversal.TraversalDescription;
import org.neo4j.graphdb.traversal.Traverser;

import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by ml on 01/04/15.
 */
public class Database {

    private static final Logger log = Logger.getLogger(Database.class.getName());

    private final static Label sampleLabel = DynamicLabel.label("Sample");
    private final static Label phenotypeLabel = DynamicLabel.label("Phenotype");
    private final static Label variantLabel = DynamicLabel.label("Variant");
    private final static Label geneLabel = DynamicLabel.label("Gene");
    private final static Label annotationLabel = DynamicLabel.label("Annotation");

    private ArrayList<Long> sampleIDNodeIDs = new ArrayList<>();
    private final GraphDatabaseService graphDb;

    public Database(final GraphDatabaseService graphDb){
        this.graphDb = graphDb;
    }

    private static enum relTypes implements RelationshipType
    {
        HAS_GENOTYPE,
        IN_GENE,
        HAS_ANNOTATION
    }
    public static void createIndexAndWait(final GraphDatabaseService graphDb, final Label label, final String property){

        IndexDefinition indexDefinition;

        try ( Transaction tx = graphDb.beginTx() )
        {
            Schema schema = graphDb.schema();
            indexDefinition = schema.indexFor( label )
                    .on(property)
                    .create();
            tx.success();
        }

        try ( Transaction tx = graphDb.beginTx() )
        {
            Schema schema = graphDb.schema();
            schema.awaitIndexOnline(indexDefinition, 10, TimeUnit.SECONDS);
        }
    }
    public static void createIndex(final GraphDatabaseService graphDb, final Label label, final String property){

        try ( Transaction tx = graphDb.beginTx() )
        {
            graphDb.schema()
                    .indexFor(label)
                    .on(property)
                    .create();

            tx.success();
        }

    }
    public static void createConstraint(final GraphDatabaseService graphDb, final Label label, final String property) {

        try ( Transaction tx = graphDb.beginTx() )
        {
            graphDb.schema()
                    .constraintFor(label)
                    .assertPropertyIsUnique(property)
                    .create();

            tx.success();
        }

    }
    public static void dropIndex(final GraphDatabaseService graphDb, final Label label){

        try ( Transaction tx = graphDb.beginTx() )
        {
            for ( IndexDefinition indexDefinition : graphDb.schema()
                    .getIndexes(label) )
            {
                indexDefinition.drop();

            }

            tx.success();
        }

    }
    public static void addNode(final GraphDatabaseService graphDb, final Label label, final String field, final Object value){

        try ( Transaction tx = graphDb.beginTx() )
        {
            Node userNode = graphDb.createNode( label );
            userNode.setProperty(field, value);

            tx.success();
        }

    }
    public static boolean hasNode(final GraphDatabaseService graphDb, final Label label, final String field, final Object value){

        try ( Transaction tx = graphDb.beginTx() )
        {
            if (graphDb.findNode(label, field, value) == null){
                return false;
            } else {
                return true;
            }
        }

    }
    public static ArrayList<Node> getNodes(final GraphDatabaseService graphDb, final Label label, final String field, final Object value){

        ArrayList<Node> nodes = new ArrayList<>();

        try ( Transaction tx = graphDb.beginTx() )
        {
            try ( ResourceIterator<Node> users = graphDb.findNodes( label, field, value ) )
            {

                while ( users.hasNext() )
                {
                    nodes.add( users.next() );
                }

            }

        }

        return nodes;
    }
    public static ArrayList<Long> getNodeIDs(final GraphDatabaseService graphDb, final Label label, final String field, final Object value){

        Node node;
        ArrayList<Long> nodeIDs = new ArrayList<>();

        try ( Transaction tx = graphDb.beginTx() )
        {
            try ( ResourceIterator<Node> users = graphDb.findNodes( label, field, value ) )
            {

                while ( users.hasNext() )
                {
                    node = users.next();
                    nodeIDs.add(node.getId());
                }

            }

        }

        return nodeIDs;
    }
    private static Traverser traverserAnnotationNodes(final GraphDatabaseService graphDb, final Node variant)
    {
        TraversalDescription td = graphDb.traversalDescription()
                .breadthFirst()
                .relationships( relTypes.HAS_ANNOTATION, Direction.OUTGOING )
                .evaluator( Evaluators.excludeStartPosition() );
        return td.traverse(variant);
    }
    private static Node getNodeById(final GraphDatabaseService graphDb, final long nodeID)
    {
        return graphDb.getNodeById( nodeID )
                .getSingleRelationship( relTypes.HAS_GENOTYPE, Direction.OUTGOING )
                .getEndNode();
    }
    private Traverser findHackers(final GraphDatabaseService graphDb, final Node startNode)
    {
        TraversalDescription td = graphDb.traversalDescription()
                .breadthFirst()
                .relationships(relTypes.HAS_GENOTYPE, Direction.OUTGOING)
                .relationships(relTypes.HAS_ANNOTATION, Direction.OUTGOING)
                .evaluator(
                        Evaluators.includeWhereLastRelationshipTypeIs(relTypes.HAS_ANNOTATION));
        return td.traverse(startNode);
    }
    public static void linkSamplesToVariants(final GraphDatabaseService graphDb, final VCFFileReader vcfFileReader){

        log.log(Level.INFO, "Linking samples to variants ...");

        //get sampleID nodeIDs
        ArrayList<Long> sampleIDNodes = new ArrayList<>();

        for (String sampleID : vcfFileReader.getFileHeader().getSampleNamesInOrder()){
            sampleIDNodes.add(getNodes(graphDb, sampleLabel, "SampleID", sampleID).get(0));
        }

        Node sampleIDNode, variantNode;
        Iterator<VariantContext> it = vcfFileReader.iterator();

        //read VCF file
        while (it.hasNext()) {
            VariantContext variant = it.next();

            //skip variants failing QC
            if (variant.isFiltered() || !variant.isVariant()) continue;

            //loop over genotypes
            for (int n = 0; n < variant.getGenotypes().size(); ++n){

                //skip wildtype genotypes
                if (!variant.getGenotypes().get(n).isHomRef() && !variant.getGenotypes().get(n).isNoCall()){
                    for (Allele allele : variant.getGenotypes().get(n).getAlleles()){

                        //make link
                        try (Transaction tx = graphDb.beginTx()) {

                            Relationship relationship = sampleIDNodes.get(n).createRelationshipTo(
                                    getNodes(graphDb, variantLabel, "VariantID", variant.getContig() + ":" + variant.getStart() + "-" + variant.getEnd() + variant.getReference() + ">" + allele.getBaseString()).get(0),
                                    relTypes.HAS_GENOTYPE);

                            if (variant.getGenotypes().get(n).isHet()) relationship.setProperty("Het", true); else relationship.setProperty("Het", false);
                            if (variant.getGenotypes().get(n).isHomVar()) relationship.setProperty("HomVar", true); else relationship.setProperty("HomVar", false);
                            if (variant.getGenotypes().get(n).isMixed()) relationship.setProperty("Mixed", true); else relationship.setProperty("Mixed", false);

                            tx.success();
                        }

                    }
                }

            }

        }

    }
    public static void addVariantsAndAnnotations(final GraphDatabaseService graphDb, final VCFFileReader vcfFileReader){

        log.log(Level.INFO, "Adding variants and annotations ...");

        Node variantNode, annotationNode;
        Iterator<VariantContext> it = vcfFileReader.iterator();
        HashSet<VEPAnnotation> splitAnnotations = new HashSet<>();

        //read VCF file
        while (it.hasNext()){
            VariantContext variant = it.next();

            //skip variants failing QC
            if (variant.isFiltered() || !variant.isVariant()) continue;

            //split annotations and make unique
            try {

                //one annotation
                VEPAnnotation splitAnnotation = new VEPAnnotation((String)variant.getAttribute("CSQ"));
                splitAnnotation.parseAnnotation();

                splitAnnotations.add(splitAnnotation);

            } catch (ClassCastException e){

                //multiple annotations
                for (String annotation : (ArrayList<String>) variant.getAttribute("CSQ")){

                    //multiple annotations
                    VEPAnnotation splitAnnotation = new VEPAnnotation(annotation);
                    splitAnnotation.parseAnnotation();

                    splitAnnotations.add(splitAnnotation);

                }

            }

            //loop over alternative alleles for this record
            for (Allele allele : variant.getAlleles()){
                if (allele.isNonReference()){

                    if (!hasNode(graphDb, variantLabel, "VariantID", variant.getContig() + ":" + variant.getStart() + "-" + variant.getEnd() + variant.getReference() + ">" + allele.getBaseString())){

                        //add variant
                        try (Transaction tx = graphDb.beginTx()) {

                            variantNode = graphDb.createNode();
                            variantNode.addLabel(variantLabel);
                            variantNode.setProperty("Contig", variant.getContig());
                            variantNode.setProperty("Start", variant.getStart());
                            variantNode.setProperty("End", variant.getEnd());
                            variantNode.setProperty("Reference", variant.getReference().getBaseString());
                            variantNode.setProperty("Alternative", allele.getBaseString());

                            tx.success();
                        }

                        //loop over annotations
                        for (VEPAnnotation annotation : splitAnnotations){

                            //skip annotations without consequences
                            if (annotation.getConsequences() == null || annotation.getConsequences().size() == 0) continue;

                            //correct annotation for this allele
                            if (allele.getBaseString().equals(annotation.getAllele())){

                                //make annotation node
                                //TODO add allele frequencies
                                try (Transaction tx = graphDb.beginTx()) {

                                    annotationNode = graphDb.createNode();
                                    annotationNode.addLabel(annotationLabel);

                                    if(annotation.getExon() != null) {
                                        String[] fields = annotation.getExon().split("/");
                                        annotationNode.setProperty("Exon", fields[0]);
                                        annotationNode.setProperty("TotalExons", fields[1]);
                                    }
                                    if(annotation.getFeature() != null) annotationNode.setProperty("Feature", annotation.getFeature());
                                    if(annotation.getFeatureType() != null) annotationNode.setProperty("FeatureType", annotation.getFeatureType());
                                    if(annotation.getGene() != null) annotationNode.setProperty("Gene", annotation.getGene());
                                    if(annotation.getHgvsCoding() != null) annotationNode.setProperty("HGVSc", annotation.getHgvsCoding());
                                    if(annotation.getHgvsProtein() != null) annotationNode.setProperty("HGVSp", annotation.getHgvsProtein());
                                    if(annotation.getIntron() != null) {
                                        String[] fields = annotation.getIntron().split("/");
                                        annotationNode.setProperty("Intron", fields[0]);
                                        annotationNode.setProperty("TotalIntrons", fields[1]);
                                    }
                                    if(annotation.getPolyphen() != null) annotationNode.setProperty("Polyphen", annotation.getPolyphen());
                                    if(annotation.getSift() != null) annotationNode.setProperty("Sift", annotation.getSift());
                                    if(annotation.getStrand() != null) annotationNode.setProperty("Strand", annotation.getStrand());
                                    if(annotation.getSymbol() != null) annotationNode.setProperty("Symbol", annotation.getSymbol());

                                    tx.success();
                                }

                                //link annotation to variant
                                try (Transaction tx = graphDb.beginTx()) {

                                    Relationship relationship = variantNode.createRelationshipTo(annotationNode, relTypes.HAS_ANNOTATION);

                                    //add consequences

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

                                }

                            }


                        } //done looping over annotations

                    }

                }

            }

            splitAnnotations.clear();
        }

    }
    public void addSampleNodes(final GraphDatabaseService graphDb, final VCFFileReader vcfFileReader){
        log.log(Level.INFO, "Adding sample nodes ...");

        Node node;

        for (String sampleID : vcfFileReader.getFileHeader().getSampleNamesInOrder()){

            if (!hasNode(graphDb, sampleLabel, "SampleID", sampleID)){

                //add sampleID
                try (Transaction tx = graphDb.beginTx()) {

                    node = graphDb.createNode();
                    node.addLabel(sampleLabel);
                    node.setProperty("SampleID", sampleID);
                    sampleIDNodeIDs.add(node.getId());

                    tx.success();
                }

            } else {
                log.log(Level.INFO, sampleID + " already exists in database.");
            }

        }

    }
    public static Label getSampleLabel() {
        return sampleLabel;
    }
    public static Label getPhenotypeLabel() {
        return phenotypeLabel;
    }
    public static Label getVariantLabel() {
        return variantLabel;
    }
    public static Label getGeneLabel() {
        return geneLabel;
    }
    public static Label getAnnotationLabel() {
        return annotationLabel;
    }

}
