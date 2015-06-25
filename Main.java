package nhs.genetics.cardiff;

import htsjdk.variant.vcf.VCFFileReader;
import org.neo4j.io.fs.FileUtils;

import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

public class Main {

    private static final Logger log = Logger.getLogger(Main.class.getName());
    private static final double version = 0.1;

    public static void main(String[] args) {

        if (args.length != 2) {
            System.err.println("Usage: <AnnotatedVCFFile> <Neo4jDBPath>");
            System.exit(1);
        }

        log.log(Level.INFO, "ImportToNeo4j v" + version);
        log.log(Level.INFO, "Importing " + args[0] + " to " + args[1]);

        //create VCF file parser
        VCFFileReader vcfFileReader = new VCFFileReader(new File(args[0]), new File(args[0] + ".idx"));

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
        variantDatabase.addSampleNodes();
        variantDatabase.addVariants();
        variantDatabase.addFunctionalAnnotations();
        variantDatabase.linkSamplesToVariants();

        variantDatabase.shutdownDatabase();

    }

}
