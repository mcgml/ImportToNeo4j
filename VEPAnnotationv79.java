package nhs.genetics.cardiff;

import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by ml on 17/06/2015.
 */

public class VEPAnnotationv79 {

    private static final Logger log = Logger.getLogger(VEPAnnotationv79.class.getName());

    private String record, allele, impact,	symbol,	gene, featureType, feature, biotype, exon, intron, hgvsCoding, hgvsProtein,
    cdnaPosition, cdsPosition, proteinPosition, aminoAcids, codons, existingVariation, distance, strand, symbolSource, hgncId,
    canonical, tsl, ccds, ensp, swissprot, trembl, uniparc, sift, polyphen, domains, gMaf, afrMaf, amrMaf, asnMaf, easMaf,
    eurMaf, sasMaf, aaMaf, eaMaf, somatic, pubmed, motifName, motifPos, highInfPos, motifScoreChange;

    private HashSet<String> consequences = new HashSet<>();
    private HashSet<String> clinSigs = new HashSet<>();

    public VEPAnnotationv79(String record) {
        this.record = record;
    }

    public void parseAnnotation() {
        String[] fields = record.split("\\|");

        try {
            if (!fields[0].equals("")) {
                this.allele = fields[0];
            }

            if (!fields[1].equals("")) {
                this.gene = fields[1];
            }

            if (!fields[2].equals("")) {
                this.feature = fields[2];
            }

            if (!fields[3].equals("")) {
                this.featureType = fields[3];
            }

            if (!fields[4].equals("")) {
                for (String consequence : fields[4].split("&")){
                    this.consequences.add(consequence);
                }
            }

            if (!fields[5].equals("")) {
                this.cdnaPosition = fields[5];
            }

            if (!fields[6].equals("")) {
                this.cdsPosition = fields[6];
            }

            if (!fields[7].equals("")) {
                this.proteinPosition = fields[7];
            }

            if (!fields[8].equals("")) {
                this.aminoAcids = fields[8];
            }

            if (!fields[9].equals("")) {
                this.codons = fields[9];
            }

            if (!fields[10].equals("")) {
                this.existingVariation = fields[10];
            }

            if (!fields[11].equals("")) {
                this.aaMaf = fields[11];
            }

            if (!fields[12].equals("")) {
                this.eaMaf = fields[12];
            }

            if (!fields[13].equals("")) {
                this.exon = fields[13];
            }

            if (!fields[14].equals("")) {
                this.intron = fields[14];
            }

            if (!fields[15].equals("")) {
                this.motifName = fields[15];
            }

            if (!fields[16].equals("")) {
                this.motifPos = fields[16];
            }

            if (!fields[17].equals("")) {
                this.highInfPos = fields[17];
            }

            if (!fields[18].equals("")) {
                this.motifScoreChange = fields[18];
            }

            if (!fields[19].equals("")) {
                this.distance = fields[19];
            }

            if (!fields[20].equals("")) {
                this.strand = fields[20];
            }

            if (!fields[21].equals("")) {
                for (String clinSig : fields[21].split("&")){
                    this.clinSigs.add(clinSig);
                }
            }

            if (!fields[22].equals("")) {
                this.canonical = fields[22];
            }

            if (!fields[23].equals("")) {
                this.symbol = fields[23];
            }

            if (!fields[24].equals("")) {
                this.symbolSource = fields[24];
            }

            if (!fields[25].equals("")) {
                String[] subFields = fields[25].split("\\(");
                this.sift = subFields[0];
            }

            if (!fields[26].equals("")) {
                String[] subFields = fields[26].split("\\(");
                this.polyphen = subFields[0];
            }

            if (!fields[27].equals("")) {
                this.gMaf = fields[27];
            }

            if (!fields[28].equals("")) {
                this.biotype = fields[28];
            }

            if (!fields[29].equals("")) {
                this.ensp = fields[29];
            }

            if (!fields[30].equals("")) {
                this.domains = fields[30];
            }

            if (!fields[31].equals("")) {
                this.ccds = fields[31];
            }

            if (!fields[32].equals("")) {
                String[] subFields = fields[32].split(":");
                this.hgvsCoding = subFields[1];
            }

            if (!fields[33].equals("")) {
                String[] subFields = fields[33].split(":");

                if (subFields[1].contains("(") && subFields[1].contains(")")){
                    subFields = subFields[1].split("\\(");
                    subFields = subFields[1].split("\\)");
                    if (subFields[0].equals("p.%3D")) this.hgvsProtein = "p.="; else this.hgvsProtein = subFields[0];
                } else {
                    if (subFields[1].equals("p.%3D")) this.hgvsProtein = "p.="; else this.hgvsProtein = subFields[1];
                }
            }

            if (!fields[34].equals("")) {
                this.afrMaf = fields[34];
            }

            if (!fields[35].equals("")) {
                this.amrMaf = fields[35];
            }

            if (!fields[36].equals("")) {
                this.asnMaf = fields[36];
            }

            if (!fields[37].equals("")) {
                this.eurMaf = fields[37];
            }

            if (!fields[38].equals("")) {
                this.pubmed = fields[38];
            }

        } catch (ArrayIndexOutOfBoundsException e) {
            log.log(Level.FINE, e.getMessage());
        }

    }

    public String getAllele() {
        return allele;
    }
    public String getGene() {
        return gene;
    }
    public String getFeature() {
        return feature;
    }
    public String getFeatureType() {
        return featureType;
    }
    public String getCdnaPosition() {
        return cdnaPosition;
    }
    public String getCdsPosition() {
        return cdsPosition;
    }
    public String getProteinPosition() {
        return proteinPosition;
    }
    public String getAminoAcids() {
        return aminoAcids;
    }
    public String getCodons() {
        return codons;
    }
    public String getExistingVariation() {
        return existingVariation;
    }
    public String getAaMaf() {
        return aaMaf;
    }
    public String getEaMaf() {
        return eaMaf;
    }
    public String getExon() {
        return exon;
    }
    public String getIntron() {
        return intron;
    }
    public String getMotifName() {
        return motifName;
    }
    public String getMotifPos() {
        return motifPos;
    }
    public String getHighInfPos() {
        return highInfPos;
    }
    public String getMotifScoreChange() {
        return motifScoreChange;
    }
    public String getDistance() {
        return distance;
    }
    public String getStrand() {
        return strand;
    }
    public String getCanonical() {
        return canonical;
    }
    public String getSymbol() {
        return symbol;
    }
    public String getSymbolSource() {
        return symbolSource;
    }
    public String getSift() {
        return sift;
    }
    public String getPolyphen() {
        return polyphen;
    }
    public String getgMaf() {
        return gMaf;
    }
    public String getBiotype() {
        return biotype;
    }
    public String getEnsp() {
        return ensp;
    }
    public String getDomains() {
        return domains;
    }
    public String getCcds() {
        return ccds;
    }
    public String getHgvsCoding() {
        return hgvsCoding;
    }
    public String getHgvsProtein() {
        return hgvsProtein;
    }
    public String getAfrMaf() {
        return afrMaf;
    }
    public String getAmrMaf() {
        return amrMaf;
    }
    public String getAsnMaf() {
        return asnMaf;
    }
    public String getEurMaf() {
        return eurMaf;
    }
    public String getPubmed() {
        return pubmed;
    }
    public HashSet<String> getConsequences() {
        return consequences;
    }
    public HashSet<String> getClinSigs() {
        return clinSigs;
    }

}