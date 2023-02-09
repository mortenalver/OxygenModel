package fishmodel.pellets;

/**
 * Super individual representation of salmon (Alver et al., 2004) with weight, numbers and stomach contents.
 */
public class SimpleFish {

    // Precomputed grouping parameters for 7 groups with boundaries at +-0.333 STD, +-1.0 STD and +-2.0 STD:
    final int nGroups = 7;
    final double[] Ndist = new double[] {0.0228, 0.1359, 0.2120, 0.2586, 0.2120, 0.1359, 0.0228};
    final double[] wDev = new double[] {-2.3577, -1.374, -0.6341, 0, 0.6341, 1.374, 2.3577};

    double a_1 = 5.2591e-6, a_2 = 0.7639; // Stomach evacuation parameters.

    double[] N, weight, V, ingested, ingRate;

    public SimpleFish(double Nfish, double meanWeight, double stdevWeight) {
        N = new double[nGroups];
        weight = new double[nGroups];
        V = new double[nGroups];
        ingested = new double[nGroups];
        ingRate = new double[nGroups];
        for (int i=0; i<nGroups; i++) {
            N[i] = Nfish*Ndist[i];
            weight[i] = meanWeight + wDev[i]*stdevWeight;
        }
    }

    public int getNGroups() {
        return nGroups;
    }

    public double getTotalN() {
        double sum = 0;
        for (int i = 0; i < N.length; i++) {
            sum += N[i];
        }
        return sum;
    }

    public double getTotalW() {
        double sum = 0;
        for (int i = 0; i < weight.length; i++) {
            sum += N[i]*weight[i];
        }
        return sum;
    }

    public double getMaxW() {
        double maxW = 0;
        for (int i = 0; i < weight.length; i++) {
            if (weight[i] > maxW)
                maxW = weight[i];
        }
        return maxW;
    }

    public double getN(int group) {
        return N[group];
    }

    public double getW(int group) {
        return weight[group];
    }

    public double getV(int group) {
        return V[group];
    }

    public void resetAllV() {
        for (int i = 0; i < V.length; i++) {
            V[i] = 0.;
        }
    }

    public double getIngested(int group) {
        return ingested[group];
    }

    public double getIngRate(int group) {
        return ingRate[group];
    }

    public void setN(int group, double N) {
        this.N[group] = N;
    }

    public void setW(int group, double W) {
        weight[group] = W;
    }

    public void addIngestion(int group, double ingestion) {
        this.V[group] += ingestion;
        this.ingested[group] += ingestion;
    }

    public void stepGutContent(int group, double dt, double T_w) {
        V[group] = V[group]*(1 - dt*a_1*Math.pow(T_w, a_2));
    }

    public void setIngRate(int group, double rate) {
        ingRate[group] = rate;
    }

    public double getAppetite(int group) {
        double mgv = getMaxGutVolume(weight[group]);
        double rsv = V[group]/mgv;
        if (rsv > 0.3)
            return Math.max(0, 0.5 - 0.57*(rsv-0.3)/(rsv-0.2));
        else
            return Math.min(1, 0.5 + 0.67*(0.3-rsv)/(0.4-rsv));
    }

    public double getMaxGutVolume(double weight) {
        // Burley and Vigg (1989):
        return 0.0007*Math.pow(weight, 1.3796);
    }
}
