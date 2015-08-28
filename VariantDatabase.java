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

    private File neo4jDBPath;
    private GraphDatabaseService graphDb;
    private VCFFileReader variantVcfFileReader, oneKgP3VcfFileReader;
    private ArrayList<VariantContext> variants = new ArrayList<>();
    private HashMap<String, HashSet<VEPAnnotation>> vepAnnotations = new HashMap<>();
    private HashMap<GenomeVariant, Node> variantNodes = new HashMap<>();
    private HashMap<String, Node> patientNodes = new HashMap<>();
    private HashMap<String, Node> sampleNodes = new HashMap<>();
    private HashMap<GenomeVariant, HashMap<String, Node>> annotationNodes = new HashMap<>();
    private HashMap<String, Node> symbolNodes = new HashMap<>();
    private HashMap<String, Node> featureNodes = new HashMap<>();
    private String libraryId = "150716_D00501_0047_BHB092ADXX";

    private Label patientLabel = DynamicLabel.label("Patient");
    private Label sampleLabel = DynamicLabel.label("Sample");
    private Label snpLabel = DynamicLabel.label("Snp");
    private Label indelLabel = DynamicLabel.label("Indel");
    private Label annotationLabel = DynamicLabel.label("Annotation");
    private Label symbolLabel = DynamicLabel.label("Symbol");
    private Label canonicalLabel = DynamicLabel.label("Canonical");
    private Label featureLabel = DynamicLabel.label("Feature");

    public VariantDatabase(VCFFileReader variantVcfFileReader, VCFFileReader oneKgP3VcfFileReader, File neo4jDBPath){
        this.variantVcfFileReader = variantVcfFileReader;
        this.neo4jDBPath = neo4jDBPath;
        this.oneKgP3VcfFileReader = oneKgP3VcfFileReader;
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
    public void startDatabase(){
        log.log(Level.INFO, "Starting database ...");

        //create DB
        graphDb = new GraphDatabaseFactory()
                .newEmbeddedDatabaseBuilder(neo4jDBPath)
                .newGraphDatabase();

        Neo4j.registerShutdownHook(graphDb);
    }
    public void createIndexes() {
        Neo4j.createConstraint(graphDb, patientLabel, "PatientId");
        Neo4j.createConstraint(graphDb, sampleLabel, "SampleId");
        Neo4j.createConstraint(graphDb, snpLabel, "VariantId");
        Neo4j.createConstraint(graphDb, indelLabel, "VariantId");
        Neo4j.createConstraint(graphDb, symbolLabel, "SymbolId");
        Neo4j.createConstraint(graphDb, featureLabel, "FeatureId");
        Neo4j.createConstraint(graphDb, canonicalLabel, "FeatureId");
    }

    public void loadVCFFile(){
        log.log(Level.INFO, "Loading VCF into memory ...");

        int n = 0;

        String variantLookup;
        Iterator<VariantContext> variantContextIterator = variantVcfFileReader.iterator();

        //read VCF file
        while (variantContextIterator.hasNext()) {
            VariantContext variant = variantContextIterator.next();

            if (variant.isFiltered() || !variant.isVariant()){
                continue;
            }

            if (n > 2000) break;
            ++n;

            variants.add(variant);

            //split annotations and make unique
            try {

                //one annotation
                VEPAnnotation vepAnnotation = new VEPAnnotation((String) variant.getAttribute("CSQ"));
                vepAnnotation.parseAnnotation();

                variantLookup = variant.getContig() + ":" + variant.getStart() + vepAnnotation.getAllele();

                //transcript, ensg, hgnc
                if (vepAnnotation.getFeature() != null && vepAnnotation.getGene() != null && vepAnnotation.getSymbol() != null) {

                    if(!vepAnnotations.containsKey(variantLookup)){
                        vepAnnotations.put(variantLookup, new HashSet<VEPAnnotation>());
                    }

                    vepAnnotations.get(variantLookup).add(vepAnnotation);

                }

            } catch (ClassCastException e) {

                //multiple annotations
                for (String annotation : (ArrayList<String>) variant.getAttribute("CSQ")) {

                    VEPAnnotation vepAnnotation = new VEPAnnotation(annotation);
                    vepAnnotation.parseAnnotation();

                    variantLookup = variant.getContig() + ":" + variant.getStart() + vepAnnotation.getAllele();

                    //transcript, ensg, hgnc
                    if (vepAnnotation.getFeature() != null && vepAnnotation.getGene() != null && vepAnnotation.getSymbol() != null) {

                        if(!vepAnnotations.containsKey(variantLookup)){
                            vepAnnotations.put(variantLookup, new HashSet<VEPAnnotation>());
                        }

                        vepAnnotations.get(variantLookup).add(vepAnnotation);

                    }

                }

            }

        }
    }

    public void addPatientNodes(){
        log.log(Level.INFO, "Adding patient nodes ...");
        //TODO
    }
    public void addSampleNodes() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Adding sample nodes ...");

        for (String sampleId : variantVcfFileReader.getFileHeader().getSampleNamesInOrder()){
            sampleNodes.put(sampleId, Neo4j.matchOrCreateUniqueNode(graphDb, sampleLabel, "SampleId", sampleId));
        }

    }
    public void addVariantNodesAndGenotypeRelationships() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Adding variants and genotypes ...");

        GenomeVariant genomeVariant;

        //read VCF records
        for (VariantContext variant : variants) {
            Iterator<Genotype> genotypeIterator = variant.getGenotypes().iterator();

            if (!variant.isVariant()) continue;

            while (genotypeIterator.hasNext()) {
                Genotype genotype = genotypeIterator.next();

                if (genotype.isNoCall() || genotype.isHomRef()) continue;
                if (genotype.getPloidy() != 2 || genotype.getAlleles().size() != 2) throw new InvalidPropertiesFormatException("Allele " + genotype.getAlleles().toString() + " is not diploid");
                if (genotype.isMixed()) continue;

                if (genotype.isHomVar()) {

                    genomeVariant = new GenomeVariant(variant.getContig(), variant.getStart(), variant.getReference().getBaseString(), genotype.getAlleles().get(1).getBaseString());
                    genomeVariant.convertToMinimalRepresentation();

                    addVariantNodesAndGenotypeRelationshipsHelper(sampleNodes.get(genotype.getSampleName()), genomeVariant, genotype.getGQ(), relTypes.HAS_HOM_VARIANT, variant.getID());

                } else if (genotype.isHetNonRef()) {

                    //het variant1
                    genomeVariant = new GenomeVariant(variant.getContig(), variant.getStart(), variant.getReference().getBaseString(), genotype.getAlleles().get(0).getBaseString());
                    genomeVariant.convertToMinimalRepresentation();

                    addVariantNodesAndGenotypeRelationshipsHelper(sampleNodes.get(genotype.getSampleName()), genomeVariant, genotype.getGQ(), relTypes.HAS_HET_VARIANT, variant.getID());

                    //het variant2
                    genomeVariant = new GenomeVariant(variant.getContig(), variant.getStart(), variant.getReference().getBaseString(), genotype.getAlleles().get(1).getBaseString());
                    genomeVariant.convertToMinimalRepresentation();

                    addVariantNodesAndGenotypeRelationshipsHelper(sampleNodes.get(genotype.getSampleName()), genomeVariant, genotype.getGQ(), relTypes.HAS_HET_VARIANT, variant.getID());

                } else if (genotype.isHet()) {

                    genomeVariant = new GenomeVariant(variant.getContig(), variant.getStart(), variant.getReference().getBaseString(), genotype.getAlleles().get(1).getBaseString());
                    genomeVariant.convertToMinimalRepresentation();

                    addVariantNodesAndGenotypeRelationshipsHelper(sampleNodes.get(genotype.getSampleName()), genomeVariant, genotype.getGQ(), relTypes.HAS_HET_VARIANT, variant.getID());

                } else {
                    throw new InvalidPropertiesFormatException("Zygosity unknown: " + variant.toString());
                }

            }

        }

    }
    public void addSymbolNodes() throws InvalidPropertiesFormatException {
        log.log(Level.INFO, "Adding symbols ...");

        for (Map.Entry<String, HashSet<VEPAnnotation>> variant : vepAnnotations.entrySet()){
            for (VEPAnnotation annotation : variant.getValue()){

                //skip symbols already imported during this session
                if (!symbolNodes.containsKey(annotation.getSymbol())){
                    symbolNodes.put(annotation.getSymbol(), Neo4j.matchOrCreateUniqueNode(graphDb, symbolLabel, "SymbolId", annotation.getSymbol()));
                }

            }
        }

    }
    public void addFeatureNodes() {
        log.log(Level.INFO, "Adding features ...");

        for (Map.Entry<String, HashSet<VEPAnnotation>> variant : vepAnnotations.entrySet()){
            for (VEPAnnotation annotation : variant.getValue()){

                //skip features already imported during this session
                if (!featureNodes.containsKey(annotation.getFeature())){

                    ArrayList<Node> nodes = Neo4j.getNodes(graphDb, featureLabel, "FeatureId", annotation.getFeature());

                    if (nodes.size() == 0) {

                        try (Transaction tx = graphDb.beginTx()) {

                            Node node = graphDb.createNode();
                            node.addLabel(featureLabel);

                            node.setProperty("FeatureId", annotation.getFeature());
                            if (annotation.getFeatureType() != null) node.setProperty("FeatureType", annotation.getFeatureType());
                            if (annotation.getBiotype() != null) node.setProperty("Biotype", annotation.getBiotype());
                            if (annotation.getCanonical() != null && annotation.getCanonical().equals("YES")) node.addLabel(canonicalLabel);
                            if (annotation.getStrand() != null) node.setProperty("Strand", annotation.getStrand());

                            if(annotation.getExon() != null) {
                                String[] fields = annotation.getExon().split("/");
                                node.setProperty("TotalExons", fields[1]);
                            }

                            if(annotation.getIntron() != null) {
                                String[] fields = annotation.getIntron().split("/");
                                node.setProperty("TotalIntrons", fields[1]);
                            }

                            featureNodes.put(annotation.getFeature(), node);

                            tx.success();
                        }

                    } else {
                        featureNodes.put(annotation.getFeature(), nodes.get(0));
                    }

                }

            }
        }

    }
    public void addFunctionalAnnotationNodes() {
        log.log(Level.INFO, "Adding functional annotations ...");

        //TODO add more annotations

        //loop over variants
        for (Map.Entry<GenomeVariant, Node> variant : variantNodes.entrySet()){
            if (vepAnnotations.containsKey(variant.getKey().getVEPVariant())){

                if (Neo4j.findNeighbourNodes(graphDb, variant.getValue(), annotationLabel, Direction.OUTGOING).size() == 0){

                    //initialise hash
                    if (!annotationNodes.containsKey(variant.getKey())) {
                        annotationNodes.put(variant.getKey(), new HashMap<String, Node>());
                    }

                    //loop over functional annotations for this variant
                    for (VEPAnnotation annotation : vepAnnotations.get(variant.getKey().getVEPVariant())) {

                        try (Transaction tx = graphDb.beginTx()) {

                            Node node = graphDb.createNode();
                            node.addLabel(annotationLabel);

                            if(annotation.getHgvsCoding() != null) node.setProperty("HGVSc", annotation.getHgvsCoding());
                            if(annotation.getHgvsProtein() != null) node.setProperty("HGVSp", annotation.getHgvsProtein());
                            if(annotation.getPolyphen() != null) node.setProperty("PolyPhen", annotation.getPolyphen());
                            if(annotation.getSift() != null) node.setProperty("SIFT", annotation.getSift());

                            if(annotation.getExon() != null) {
                                String[] fields = annotation.getExon().split("/");
                                node.setProperty("Exon", fields[0]);
                            }

                            if(annotation.getIntron() != null) {
                                String[] fields = annotation.getIntron().split("/");
                                node.setProperty("Intron", fields[0]);
                            }

                            annotationNodes.get(variant.getKey()).put(annotation.getFeature(), node);

                            tx.success();
                        }

                    }
                }

            } else {
                log.log(Level.WARNING, "No annotation for: " + variant.getKey().getConcatenatedVariant());
            }
        }

    }
    public void addConsequenceRelationships() {
        log.log(Level.INFO, "Linking variants To functional annotations ...");

        //loop over variants
        for (Map.Entry<GenomeVariant, Node> variant : variantNodes.entrySet()){
            if (vepAnnotations.containsKey(variant.getKey().getVEPVariant())){

                HashMap<String, Object> properties = new HashMap<>();

                //loop over all annotations by for this allele
                for (VEPAnnotation annotation : vepAnnotations.get(variant.getKey().getVEPVariant())){

                    if (annotation.getConsequences().contains("transcript_ablation")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_TRANSCRIPT_ABLATION_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("splice_acceptor_variant")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_SPLICE_ACCEPTOR_VARIANT_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("splice_donor_variant")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_SPLICE_DONOR_VARIANT_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("stop_gained")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_STOP_GAINED_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("frameshift_variant")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_FRAMESHIFT_VARIANT_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("stop_lost")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_STOP_LOST_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("start_lost")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_START_LOST_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("transcript_amplification")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_TRANSCRIPT_AMPLIFICATION_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("inframe_insertion")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_INFRAME_INSERTION_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("inframe_deletion")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_INFRAME_DELETION_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("missense_variant")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_MISSENSE_VARIANT_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("protein_altering_variant")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_PROTEIN_ALTERING_VARIANT_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("splice_region_variant")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_SPLICE_REGION_VARIANT_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("incomplete_terminal_codon_variant")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_INCOMPLETE_TERMINAL_CODON_VARIANT_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("stop_retained_variant")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_STOP_RETAINED_VARIANT_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("synonymous_variant")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_SYNONYMOUS_VARIANT_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("coding_sequence_variant")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_CODING_SEQUENCE_VARIANT_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("mature_miRNA_variant")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_MATURE_MIRNA_VARIANT_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("5_prime_UTR_variant")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_5_PRIME_UTR_VARIANT_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("3_prime_UTR_variant")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_3_PRIME_UTR_VARIANT_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("non_coding_transcript_exon_variant")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_NON_CODING_TRANSCRIPT_EXON_VARIANT_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("intron_variant")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_INTRON_VARIANT_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("NMD_transcript_variant")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_NMD_TRANSCRIPT_VARIANT_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("non_coding_transcript_variant")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_NON_CODING_TRANSCRIPT_VARIANT_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("upstream_gene_variant")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_UPSTREAM_GENE_VARIANT_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("downstream_gene_variant")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_DOWNSTREAM_GENE_VARIANT_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("TFBS_ablation")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_TFBS_ABLATION_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("TFBS_amplification")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_TFBS_AMPLIFICATION_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("TF_binding_site_variant")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_TF_BINDING_SITE_VARIANT_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("regulatory_region_ablation")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_REGULATORY_REGION_ABLATION_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("regulatory_region_amplification")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_REGULATORY_REGION_AMPLIFICATION_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("feature_elongation")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_FEATURE_ELONGATION_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("regulatory_region_variant")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_REGULATORY_REGION_VARIANT_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("feature_truncation")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_FEATURE_TRUNCATION_CONSEQUENCE, properties, false);
                    if (annotation.getConsequences().contains("intergenic_variant")) Neo4j.createRelationship(graphDb, variantNodes.get(variant.getKey()), annotationNodes.get(variant.getKey()).get(annotation.getFeature()), relTypes.HAS_INTERGENIC_VARIANT_CONSEQUENCE, properties, false);

                }
            }

        }
    }
    public void addInFeatureRelationships() {
        log.log(Level.INFO, "Linking annotations to features ...");

        //loop over variants
        for (Map.Entry<GenomeVariant, Node> variant : variantNodes.entrySet()){
            if (vepAnnotations.containsKey(variant.getKey().getVEPVariant())){

                //loop over functional annotations for this variant
                for (VEPAnnotation annotation : vepAnnotations.get(variant.getKey().getVEPVariant())) {
                    Neo4j.createRelationship(graphDb, annotationNodes.get(variant.getKey()).get(annotation.getFeature()), featureNodes.get(annotation.getFeature()), relTypes.IN_FEATURE, new HashMap<String, Object>(), false);
                }

            }
        }

    }
    public void addInSymbolRelationships() {
        log.log(Level.INFO, "Linking features to symbol ...");

        //loop over variants
        for (Map.Entry<GenomeVariant, Node> variant : variantNodes.entrySet()) {

            if (vepAnnotations.containsKey(variant.getKey().getVEPVariant())) {

                //loop over functional annotations for this variant
                for (VEPAnnotation annotation : vepAnnotations.get(variant.getKey().getVEPVariant())) {
                    Neo4j.createRelationship(graphDb, featureNodes.get(annotation.getFeature()), symbolNodes.get(annotation.getSymbol()), relTypes.IN_SYMBOL, new HashMap<String, Object>(), false);
                }
            }

        }

    }

    private void addVariantNodesAndGenotypeRelationshipsHelper(Node sampleNode, GenomeVariant genomeVariant, int genotypeQuality, RelationshipType relationshipType, String rsId) throws InvalidPropertiesFormatException {

        boolean hasGenotype = false;

        //add variant node if not already during this session
        if (!variantNodes.containsKey(genomeVariant)){

            if (genomeVariant.getRef().length() == 1 && genomeVariant.getAlt().length() == 1) {

                ArrayList<Node> nodes = Neo4j.getNodes(graphDb, snpLabel, "VariantId", genomeVariant.getConcatenatedVariant());

                //variant not present in DB -- add variant
                if (nodes.size() == 0){

                    HashMap<String, Object> properties = getPopulationFrequencies(genomeVariant);
                    properties.put("VariantId", genomeVariant.getConcatenatedVariant());
                    properties.put("rsId", rsId);

                    variantNodes.put(genomeVariant, Neo4j.addNode(graphDb, snpLabel, properties));
                }

            } else if (genomeVariant.getRef().length() != 1 || genomeVariant.getAlt().length() != 1) {

                ArrayList<Node> nodes = Neo4j.getNodes(graphDb, indelLabel, "VariantId", genomeVariant.getConcatenatedVariant());

                //variant not present in DB -- add variant
                if (nodes.size() == 0){

                    HashMap<String, Object> properties = getPopulationFrequencies(genomeVariant);
                    properties.put("VariantId", genomeVariant.getConcatenatedVariant());
                    properties.put("rsId", rsId);

                    variantNodes.put(genomeVariant, Neo4j.addNode(graphDb, indelLabel, properties));
                }

            } else {
                throw new InvalidPropertiesFormatException("Variant " + genomeVariant.getConcatenatedVariant() + " not SNP or Indel");
            }

        }

        //check if this genotype already exists in the DB
        try ( Transaction tx = graphDb.beginTx() ){
            for (Relationship relationship : sampleNode.getRelationships(Direction.OUTGOING)){

                if (relationship.getOtherNode(sampleNode) == variantNodes.get(genomeVariant) && relationship.getProperty("LibraryId").equals(libraryId)){
                    hasGenotype = true;
                    break;
                }

            }
        }

        if (!hasGenotype){
            HashMap<String, Object> properties = new HashMap<>();

            properties.put("GQ", genotypeQuality);
            properties.put("LibraryId", libraryId);

            Neo4j.createRelationship(graphDb, sampleNode, variantNodes.get(genomeVariant), relationshipType, properties, true);
        }

    }
    private HashMap<String, Object> getPopulationFrequencies(GenomeVariant genomeVariant){

        //TODO check vars are lookup correctly

        GenomeVariant tempVariant;
        HashMap<String, Object> populationFrequencies = new HashMap<>();
        ArrayList<String> infoAtrributes;

        Iterator<VariantContext> variantContextIterator = oneKgP3VcfFileReader.query(genomeVariant.getContig(), genomeVariant.getPos(), genomeVariant.getPos());

        //get 1kg variants
        while (variantContextIterator.hasNext()) {
            VariantContext vcfVariant = variantContextIterator.next();

            for (int n = 0; n < vcfVariant.getAlternateAlleles().size(); ++n){

                tempVariant = new GenomeVariant(vcfVariant.getContig(), vcfVariant.getStart(), vcfVariant.getReference().getBaseString(), vcfVariant.getAlternateAllele(n).getBaseString());
                tempVariant.convertToMinimalRepresentation();

                if (tempVariant.equals(genomeVariant)){

                    try{

                        //found variant
                        populationFrequencies.put("EAS_AF", Double.parseDouble((String) vcfVariant.getAttribute("EAS_AF")));
                        populationFrequencies.put("EUR_AF", Double.parseDouble((String) vcfVariant.getAttribute("EUR_AF")));
                        populationFrequencies.put("AFR_AF", Double.parseDouble((String) vcfVariant.getAttribute("AFR_AF")));
                        populationFrequencies.put("AMR_AF", Double.parseDouble((String) vcfVariant.getAttribute("AMR_AF")));
                        populationFrequencies.put("SAS_AF", Double.parseDouble((String) vcfVariant.getAttribute("SAS_AF")));

                    } catch (ClassCastException e) {

                        //found variant
                        infoAtrributes = (ArrayList<String>) vcfVariant.getAttribute("EAS_AF");
                        populationFrequencies.put("EAS_AF", Double.parseDouble((String) infoAtrributes.get(n)));

                        infoAtrributes = (ArrayList<String>) vcfVariant.getAttribute("EUR_AF");
                        populationFrequencies.put("EUR_AF", Double.parseDouble((String) infoAtrributes.get(n)));

                        infoAtrributes = (ArrayList<String>) vcfVariant.getAttribute("AFR_AF");
                        populationFrequencies.put("AFR_AF", Double.parseDouble((String) infoAtrributes.get(n)));

                        infoAtrributes = (ArrayList<String>) vcfVariant.getAttribute("AMR_AF");
                        populationFrequencies.put("AMR_AF", Double.parseDouble((String) infoAtrributes.get(n)));

                        infoAtrributes = (ArrayList<String>) vcfVariant.getAttribute("SAS_AF");
                        populationFrequencies.put("SAS_AF", Double.parseDouble((String) infoAtrributes.get(n)));

                    }

                    return populationFrequencies;
                }

            }

        }

        log.log(Level.WARNING, genomeVariant.getConcatenatedVariant() + " not found in poly VCF");

        return populationFrequencies;
    }

    public void shutdownDatabase(){
        log.log(Level.INFO, "Shutting down database ...");
        graphDb.shutdown();
    }

}