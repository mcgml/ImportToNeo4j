package nhs.genetics.cardiff;

import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

public class Main {

    private static final Logger log = Logger.getLogger(Main.class.getName());

    public static void main(String[] args) {

        if (args.length != 2) {
            System.err.println("Usage: <AnnotatedVCFFile> <Neo4jDBPath>");
            System.exit(1);
        }

        //create VCF file parser
        VCFFileReader vcfFileReader = new VCFFileReader(new File(args[0]), new File(args[0] + ".idx"));

        //create database object
        VariantDatabase variantDatabase = new VariantDatabase(vcfFileReader, new File(args[1]));

        //delete existing DB
        try{
            variantDatabase.deleteDatabase();
        } catch (IOException e){
            log.log(Level.SEVERE, "Could not delete database: " + e.getMessage());
        }

        //create new DB
        variantDatabase.createDatabase();
        variantDatabase.addSampleNodes();




    }

}
