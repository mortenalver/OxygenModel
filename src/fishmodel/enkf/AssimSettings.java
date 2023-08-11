package fishmodel.enkf;

public class AssimSettings {

    public boolean useTwin = false; // If running EnKF, this variable is set to true if running a twin experiment
    // where measurements are acquired from the last (N-1) parallel model
    public boolean perturbTwin = false;
    public boolean perturbThisMember = false; // Perturbations will be activated if we are running with EnKF and
    // unless this member is the twin model in a twin experiment setup.
    public boolean dryRun = false; // If set to true, EnKF corrections will be computed, but not applied

    public boolean usePerturbations = !dryRun; // If false, no perturbations. IF true, perturbations will be activated where
    // appropriate (when running in MPI mode and if this process is not the twin model).

    // Use the following variables to temporarily disable EnKF corrections within a given time period.
    // There is a method further down that checks whether we are in the dropout interval for a certain time.
    public boolean correctionsDropOut = true;
    public double correctionsDropOutStartS = 8*3600,
        correctionsDropOutEndS = 10*3600;

    public int enKFInterval = 60;
    public double locDist=15/*25*/, locZMultiplier=3;

    /* Ensemble inflation "blows up" the variability within the ensemble after each analysis step. It is a simple
     * way of increasing the variances/covariances with the effect of making KF corrections stronger.
     * The boolean variable turns it on or off, and the inflation factor controls how much deviations from the
     * ensemble mean are multiplied at each analysis step.
    */
    public boolean ensembleInflation = false;
    public double ensembleInflationFactor = 1.05;
    public double ambientO2Std = 2.*0.25, ambientO2Beta = 0.2*0.05;

    public double currentStd = 0.03, currentBeta = 0.2*0.05;
    public double o2ConsStd = 0*0.25, o2ConsBeta = 0.2*0.05;

    // Parameter estimation:
    public int nPar = 4;
    double parStd0 = 0.002;
    public double[] parStd = new double[] {parStd0, parStd0, parStd0, 0.25*parStd0};

    public double[] parDepths = new double[] {5, 10, 15};


    public boolean isDropOutActive(double timeS) {
        return correctionsDropOut && (timeS>correctionsDropOutStartS) && (timeS< correctionsDropOutEndS);
    }
}
