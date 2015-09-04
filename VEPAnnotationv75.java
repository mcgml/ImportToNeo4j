package nhs.genetics.cardiff;

import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by ml on 17/06/2015.
 */
public class VEPAnnotationv75 {

    private static final Logger log = Logger.getLogger(Main.class.getName());

    private String record, allele, gene, feature, featureType, cdnaPosition, cdsPosition, proteinPosition, aminoAcids, codons,
            existingVariation, aaMaf, eaMaf, exon, intron, motifName, motifPos, highInfPos, motifScoreChange, distance, strand,
            canonical, symbol, symbolSource, sift, polyphen, gMaf, biotype, ensp, domains, ccds, hgvsCoding, hgvsProtein, afrMaf,
            amrMaf, asnMaf, eurMaf, pubmed;
    private HashSet<String> consequences = new HashSet<>();
    private HashSet<String> clinSigs = new HashSet<>();

    public VEPAnnotationv75(String record) {
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

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        VEPAnnotationv75 that = (VEPAnnotationv75) o;

        if (allele != null ? !allele.equals(that.allele) : that.allele != null) return false;
        if (gene != null ? !gene.equals(that.gene) : that.gene != null) return false;
        if (feature != null ? !feature.equals(that.feature) : that.feature != null) return false;
        if (featureType != null ? !featureType.equals(that.featureType) : that.featureType != null) return false;
        if (cdnaPosition != null ? !cdnaPosition.equals(that.cdnaPosition) : that.cdnaPosition != null) return false;
        if (cdsPosition != null ? !cdsPosition.equals(that.cdsPosition) : that.cdsPosition != null) return false;
        if (proteinPosition != null ? !proteinPosition.equals(that.proteinPosition) : that.proteinPosition != null)
            return false;
        if (aminoAcids != null ? !aminoAcids.equals(that.aminoAcids) : that.aminoAcids != null) return false;
        if (codons != null ? !codons.equals(that.codons) : that.codons != null) return false;
        if (existingVariation != null ? !existingVariation.equals(that.existingVariation) : that.existingVariation != null)
            return false;
        if (aaMaf != null ? !aaMaf.equals(that.aaMaf) : that.aaMaf != null) return false;
        if (eaMaf != null ? !eaMaf.equals(that.eaMaf) : that.eaMaf != null) return false;
        if (exon != null ? !exon.equals(that.exon) : that.exon != null) return false;
        if (intron != null ? !intron.equals(that.intron) : that.intron != null) return false;
        if (motifName != null ? !motifName.equals(that.motifName) : that.motifName != null) return false;
        if (motifPos != null ? !motifPos.equals(that.motifPos) : that.motifPos != null) return false;
        if (highInfPos != null ? !highInfPos.equals(that.highInfPos) : that.highInfPos != null) return false;
        if (motifScoreChange != null ? !motifScoreChange.equals(that.motifScoreChange) : that.motifScoreChange != null)
            return false;
        if (distance != null ? !distance.equals(that.distance) : that.distance != null) return false;
        if (strand != null ? !strand.equals(that.strand) : that.strand != null) return false;
        if (canonical != null ? !canonical.equals(that.canonical) : that.canonical != null) return false;
        if (symbol != null ? !symbol.equals(that.symbol) : that.symbol != null) return false;
        if (symbolSource != null ? !symbolSource.equals(that.symbolSource) : that.symbolSource != null) return false;
        if (sift != null ? !sift.equals(that.sift) : that.sift != null) return false;
        if (polyphen != null ? !polyphen.equals(that.polyphen) : that.polyphen != null) return false;
        if (gMaf != null ? !gMaf.equals(that.gMaf) : that.gMaf != null) return false;
        if (biotype != null ? !biotype.equals(that.biotype) : that.biotype != null) return false;
        if (ensp != null ? !ensp.equals(that.ensp) : that.ensp != null) return false;
        if (domains != null ? !domains.equals(that.domains) : that.domains != null) return false;
        if (ccds != null ? !ccds.equals(that.ccds) : that.ccds != null) return false;
        if (hgvsCoding != null ? !hgvsCoding.equals(that.hgvsCoding) : that.hgvsCoding != null) return false;
        if (hgvsProtein != null ? !hgvsProtein.equals(that.hgvsProtein) : that.hgvsProtein != null) return false;
        if (afrMaf != null ? !afrMaf.equals(that.afrMaf) : that.afrMaf != null) return false;
        if (amrMaf != null ? !amrMaf.equals(that.amrMaf) : that.amrMaf != null) return false;
        if (asnMaf != null ? !asnMaf.equals(that.asnMaf) : that.asnMaf != null) return false;
        if (eurMaf != null ? !eurMaf.equals(that.eurMaf) : that.eurMaf != null) return false;
        if (pubmed != null ? !pubmed.equals(that.pubmed) : that.pubmed != null) return false;
        if (consequences != null ? !consequences.equals(that.consequences) : that.consequences != null) return false;
        return !(clinSigs != null ? !clinSigs.equals(that.clinSigs) : that.clinSigs != null);

    }

    @Override
    public int hashCode() {
        int result = allele != null ? allele.hashCode() : 0;
        result = 31 * result + (gene != null ? gene.hashCode() : 0);
        result = 31 * result + (feature != null ? feature.hashCode() : 0);
        result = 31 * result + (featureType != null ? featureType.hashCode() : 0);
        result = 31 * result + (cdnaPosition != null ? cdnaPosition.hashCode() : 0);
        result = 31 * result + (cdsPosition != null ? cdsPosition.hashCode() : 0);
        result = 31 * result + (proteinPosition != null ? proteinPosition.hashCode() : 0);
        result = 31 * result + (aminoAcids != null ? aminoAcids.hashCode() : 0);
        result = 31 * result + (codons != null ? codons.hashCode() : 0);
        result = 31 * result + (existingVariation != null ? existingVariation.hashCode() : 0);
        result = 31 * result + (aaMaf != null ? aaMaf.hashCode() : 0);
        result = 31 * result + (eaMaf != null ? eaMaf.hashCode() : 0);
        result = 31 * result + (exon != null ? exon.hashCode() : 0);
        result = 31 * result + (intron != null ? intron.hashCode() : 0);
        result = 31 * result + (motifName != null ? motifName.hashCode() : 0);
        result = 31 * result + (motifPos != null ? motifPos.hashCode() : 0);
        result = 31 * result + (highInfPos != null ? highInfPos.hashCode() : 0);
        result = 31 * result + (motifScoreChange != null ? motifScoreChange.hashCode() : 0);
        result = 31 * result + (distance != null ? distance.hashCode() : 0);
        result = 31 * result + (strand != null ? strand.hashCode() : 0);
        result = 31 * result + (canonical != null ? canonical.hashCode() : 0);
        result = 31 * result + (symbol != null ? symbol.hashCode() : 0);
        result = 31 * result + (symbolSource != null ? symbolSource.hashCode() : 0);
        result = 31 * result + (sift != null ? sift.hashCode() : 0);
        result = 31 * result + (polyphen != null ? polyphen.hashCode() : 0);
        result = 31 * result + (gMaf != null ? gMaf.hashCode() : 0);
        result = 31 * result + (biotype != null ? biotype.hashCode() : 0);
        result = 31 * result + (ensp != null ? ensp.hashCode() : 0);
        result = 31 * result + (domains != null ? domains.hashCode() : 0);
        result = 31 * result + (ccds != null ? ccds.hashCode() : 0);
        result = 31 * result + (hgvsCoding != null ? hgvsCoding.hashCode() : 0);
        result = 31 * result + (hgvsProtein != null ? hgvsProtein.hashCode() : 0);
        result = 31 * result + (afrMaf != null ? afrMaf.hashCode() : 0);
        result = 31 * result + (amrMaf != null ? amrMaf.hashCode() : 0);
        result = 31 * result + (asnMaf != null ? asnMaf.hashCode() : 0);
        result = 31 * result + (eurMaf != null ? eurMaf.hashCode() : 0);
        result = 31 * result + (pubmed != null ? pubmed.hashCode() : 0);
        result = 31 * result + (consequences != null ? consequences.hashCode() : 0);
        result = 31 * result + (clinSigs != null ? clinSigs.hashCode() : 0);
        return result;
    }
}