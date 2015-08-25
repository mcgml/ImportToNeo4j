package nhs.genetics.cardiff;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import org.neo4j.io.fs.FileUtils;

import java.io.File;
import java.io.IOException;
import java.util.InvalidPropertiesFormatException;
import java.util.logging.Level;
import java.util.logging.Logger;

public class Main {

    private static final Logger log = Logger.getLogger(Main.class.getName());
    private static final double version = 0.1;

    public static void main(String[] args) throws InvalidPropertiesFormatException {

        if (args.length != 2) {
            System.err.println("Usage: <VEPAnnotatedVCFFile> <Neo4jDBPath>");
            System.exit(1);
        }

        log.log(Level.INFO, "ImportToNeo4j v" + version);
        log.log(Level.INFO, "Importing " + args[0] + " to " + args[1]);

        //create VCF file parser
        VCFFileReader vcfFileReader = new VCFFileReader(new File(args[0]), new File(args[0] + ".idx"));
        VCFCodec codec = new VCFCodec();

        try (AbstractFeatureReader reader = AbstractFeatureReader.getFeatureReader(args[0], codec)){
            Iterable<VCF> iter = reader.iterator();

        } catch (IOException e){
            e.printStackTrace();
        }

        //create database object
        VariantDatabase variantDatabase = new VariantDatabase(vcfFileReader, new File(args[1]));

        //delete existing DB
        try{
            FileUtils.deleteRecursively(new File(args[1]));
        } catch (IOException e){
            log.log(Level.SEVERE, "Could not delete database: " + e.getMessage());
            System.exit(1);
        }

        //create new DB
        variantDatabase.createDatabase();
        variantDatabase.createIndexes();

        variantDatabase.loadVCF();

        variantDatabase.addPatientNodes();
        variantDatabase.addSampleNodes();
        variantDatabase.addVariantsNodes();
        variantDatabase.addGenotypeRelationships();
        variantDatabase.addSymbolNodes();
        variantDatabase.addFeatureNodes();
        variantDatabase.addFunctionalAnnotationNodes();
        variantDatabase.addConsequenceRelationships();
        variantDatabase.addInFeatureRelationships();
        variantDatabase.addInSymbolRelationships();

        variantDatabase.shutdownDatabase();
        vcfFileReader.close();

    }

}
