package fishmodel.enkf;

public class AssimSettings {

    public boolean useTwin = false; // If running EnKF, this variable is set to true if running a twin experiment
    // where measurements are acquired from the last (N-1) parallel model
    public boolean perturbTwin = false;
    public boolean perturbThisMember = false; // Perturbations will be activated if we are running with EnKF and
    // unless this member is the twin model in a twin experiment setup.
    public boolean dryRun = false; // If set to true, EnKF corrections will be computed, but not applied

    public int enKFInterval = 60;
    public double locDist=25/*17*/, locZMultiplier=3;
    public double ambientO2Std = 0.1, ambientO2Beta = 0.05;

    public double currentStd = 0.02, currentBeta = 0.05;

    // Parameter estimation:
    public int nPar = 1;
    public double[] parStd = new double[] {0.002};
}
