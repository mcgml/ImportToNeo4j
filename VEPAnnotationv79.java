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
                for (String consequence : fields[1].split("&")){
                    this.consequences.add(consequence);
                }
            }
            if (!fields[2].equals("")) {
                this.impact = fields[2];
            }
            if (!fields[3].equals("")) {
                this.symbol = fields[3];
            }
            if (!fields[4].equals("")) {
                this.gene = fields[4];
            }
            if (!fields[5].equals("")) {
                this.featureType = fields[5];
            }
            if (!fields[6].equals("")) {
                this.feature = fields[6];
            }
            if (!fields[7].equals("")) {
                this.biotype = fields[7];
            }
            if (!fields[8].equals("")) {
                this.exon = fields[8];
            }
            if (!fields[9].equals("")) {
                this.intron = fields[9];
            }
            if (!fields[10].equals("")) {
                String[] subFields = fields[10].split(":");
                this.hgvsCoding = subFields[1];
            }
            if (!fields[11].equals("")) {
                String[] subFields = fields[11].split(":");

                if (subFields[1].contains("(") && subFields[1].contains(")")){
                    subFields = subFields[1].split("\\(");
                    subFields = subFields[1].split("\\)");
                    if (subFields[0].equals("p.%3D")) this.hgvsProtein = "p.="; else this.hgvsProtein = subFields[0];
                } else {
                    if (subFields[1].equals("p.%3D")) this.hgvsProtein = "p.="; else this.hgvsProtein = subFields[1];
                }
            }
            if (!fields[12].equals("")) {
                this.cdnaPosition = fields[12];
            }
            if (!fields[13].equals("")) {
                this.cdsPosition = fields[13];
            }
            if (!fields[14].equals("")) {
                this.proteinPosition = fields[14];
            }
            if (!fields[15].equals("")) {
                this.aminoAcids = fields[15];
            }
            if (!fields[16].equals("")) {
                this.codons = fields[16];
            }
            if (!fields[17].equals("")) {
                this.existingVariation = fields[17];
            }
            if (!fields[18].equals("")) {
                this.distance = fields[18];
            }
            if (!fields[19].equals("")) {
                this.strand = fields[19];
            }
            if (!fields[20].equals("")) {
                this.symbolSource = fields[20];
            }
            if (!fields[21].equals("")) {
                this.hgncId = fields[21];
            }
            if (!fields[22].equals("")) {
                this.canonical = fields[22];
            }
            if (!fields[23].equals("")) {
                this.tsl = fields[23];
            }
            if (!fields[24].equals("")) {
                this.ccds = fields[24];
            }
            if (!fields[25].equals("")) {
                this.ensp = fields[25];
            }
            if (!fields[26].equals("")) {
                this.swissprot = fields[26];
            }
            if (!fields[27].equals("")) {
                this.trembl = fields[27];
            }
            if (!fields[28].equals("")) {
                this.uniparc = fields[28];
            }
            if (!fields[29].equals("")) {
                String[] subFields = fields[29].split("\\(");
                this.sift = subFields[0];
            }
            if (!fields[30].equals("")) {
                String[] subFields = fields[30].split("\\(");
                this.polyphen = subFields[0];
            }
            if (!fields[31].equals("")) {
                this.domains = fields[31];
            }
            if (!fields[32].equals("")) {
                this.gMaf = fields[32];
            }
            if (!fields[33].equals("")) {
                this.afrMaf = fields[33];
            }
            if (!fields[34].equals("")) {
                this.amrMaf = fields[34];
            }
            if (!fields[35].equals("")) {
                this.asnMaf = fields[35];
            }
            if (!fields[36].equals("")) {
                this.easMaf = fields[36];
            }
            if (!fields[37].equals("")) {
                this.eurMaf = fields[37];
            }
            if (!fields[38].equals("")) {
                this.sasMaf = fields[38];
            }
            if (!fields[39].equals("")) {
                this.aaMaf = fields[39];
            }
            if (!fields[40].equals("")) {
                this.eaMaf = fields[40];
            }
            if (!fields[41].equals("")) {
                for (String clinSig : fields[41].split("&")){
                    this.clinSigs.add(clinSig);
                }
            }
            if (!fields[42].equals("")) {
                this.somatic = fields[42];
            }
            if (!fields[43].equals("")) {
                this.pubmed = fields[43];
            }
            if (!fields[44].equals("")) {
                this.motifName = fields[44];
            }
            if (!fields[45].equals("")) {
                this.motifPos = fields[45];
            }
            if (!fields[46].equals("")) {
                this.highInfPos = fields[46];
            }
            if (!fields[47].equals("")) {
                this.motifScoreChange = fields[47];
            }

        } catch (ArrayIndexOutOfBoundsException e) {
            log.log(Level.FINE, e.getMessage());
        }

    }

    public String getAllele() {
        return allele;
    }
    public String getImpact() {
        return impact;
    }
    public String getSymbol() {
        return symbol;
    }
    public String getGene() {
        return gene;
    }
    public String getFeatureType() {
        return featureType;
    }
    public String getFeature() {
        return feature;
    }
    public String getBiotype() {
        return biotype;
    }
    public String getExon() {
        return exon;
    }
    public String getIntron() {
        return intron;
    }
    public String getHgvsCoding() {
        return hgvsCoding;
    }
    public String getHgvsProtein() {
        return hgvsProtein;
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
    public String getDistance() {
        return distance;
    }
    public String getStrand() {
        return strand;
    }
    public String getSymbolSource() {
        return symbolSource;
    }
    public String getHgncId() {
        return hgncId;
    }
    public String getCanonical() {
        return canonical;
    }
    public String getTsl() {
        return tsl;
    }
    public String getCcds() {
        return ccds;
    }
    public String getEnsp() {
        return ensp;
    }
    public String getSwissprot() {
        return swissprot;
    }
    public String getTrembl() {
        return trembl;
    }
    public String getUniparc() {
        return uniparc;
    }
    public String getSift() {
        return sift;
    }
    public String getPolyphen() {
        return polyphen;
    }
    public String getDomains() {
        return domains;
    }
    public String getgMaf() {
        return gMaf;
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
    public String getEasMaf() {
        return easMaf;
    }
    public String getEurMaf() {
        return eurMaf;
    }
    public String getSasMaf() {
        return sasMaf;
    }
    public String getAaMaf() {
        return aaMaf;
    }
    public String getEaMaf() {
        return eaMaf;
    }
    public String getSomatic() {
        return somatic;
    }
    public String getPubmed() {
        return pubmed;
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
    public HashSet<String> getConsequences() {
        return consequences;
    }
    public HashSet<String> getClinSigs() {
        return clinSigs;
    }

}