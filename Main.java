package nhs.genetics.cardiff;

import htsjdk.variant.vcf.VCFFileReader;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import org.neo4j.io.fs.FileUtils;

import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

public class Main {

    private static final Logger log = Logger.getLogger(Main.class.getName());

    private static final File dbFilePath = new File ("/Users/ml/Documents/Neo4j/default.graphdb");

    public static void main(String[] args) throws IOException {

        if (args.length != 1) {
            System.err.println("Usage: <AnnotatedVCFFile>");
            System.exit(1);
        }

        log.log(Level.INFO, "Starting database ...");

        //delete pre-existing DB
        FileUtils.deleteRecursively(dbFilePath);

        //initalise VCF file reader
        VCFFileReader vcfFileReader = new VCFFileReader(new File(args[0]), new File(args[0] + ".idx"));

        //create DB
        GraphDatabaseService graphDb = new GraphDatabaseFactory()
                .newEmbeddedDatabaseBuilder(dbFilePath)
                .newGraphDatabase();
        registerShutdownHook(graphDb);

        //create unique indexes
        Database.createConstraint(graphDb, Database.getVariantLabel(), "VariantID");
        Database.createConstraint(graphDb, Database.getSampleLabel(), "SampleID");

        //populate database
        Database.addVariantsAndAnnotations(graphDb, vcfFileReader);
        Database.addSampleNodes(graphDb, vcfFileReader);
        Database.linkSamplesToVariants(graphDb, vcfFileReader);

        log.log(Level.INFO, "Shutting down database ...");
        graphDb.shutdown();
    }

    private static void registerShutdownHook( final GraphDatabaseService graphDb )
    {
        // Registers a shutdown hook for the Neo4j instance so that it
        // shuts down nicely when the VM exits (even if you "Ctrl-C" the
        // running application).
        Runtime.getRuntime().addShutdownHook( new Thread()
        {
            @Override
            public void run()
            {
                graphDb.shutdown();
            }
        } );
    }

}
