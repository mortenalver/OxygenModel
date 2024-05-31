package fishmodel.enkf;

public class AssimSettings {

    public boolean useEnOI = true; // Use EnOI instead of EnKF. This value is disregarded if running in MPI mode.
    public boolean hybrid_EnKF_ENOI = false;
    public double hybrid_ENOI_weight = 0.5; // Relative weighting of ENOI K matrix in hybrid setup
    public String[] enOIEnsembleFile = new String[]{
            //"C:\\Users\\alver\\Work\\BjoroyaSim\\highstorage27_06.nc"
            "static22_06.nc",
            "static23_06.nc",
            "static24_06.nc",
            "static25_06.nc",
            "static26_06.nc",
            "static27_06.nc",
            "static28_06.nc",
            "static29_06.nc",
            "static30_06.nc"
    };


    //public String enOIEnsembleFile = "C:\\Users\\alver\\Work\\BjoroyaSim\\enoitest_dr24_06.nc";
    public double enoiAlpha = 0.1;

    public boolean useTwin = false; // If running EnKF, this variable is set to true if running a twin experiment
    public String twinDataPathPrefix = "twin2";

    // where measurements are acquired from the last (N-1) parallel model
    public boolean perturbTwin = false;
    public boolean perturbThisMember = false; // Perturbations will be activated if we are running with EnKF and
    // unless this member is the twin model in a twin experiment setup.
    public boolean dryRun = false; // If set to true, EnKF corrections will be computed, but not applied

    public boolean usePerturbations = !dryRun; // If false, no perturbations. IF true, perturbations will be activated where
    // appropriate (when running in MPI mode and if this process is not the twin model).

    // Use the following variables to temporarily disable EnKF corrections within a given time period.
    // There is a method further down that checks whether we are in the dropout interval for a certain time.
    public boolean correctionsDropOut = false;
    public double correctionsDropOutStartS = 1,//8*3600,
        correctionsDropOutEndS = 6*3600;

    public int assimInterval = 60;
    public double locDist= 30/*15*/ /*30*/, locZMultiplier=3;

    /* Ensemble inflation "blows up" the variability within the ensemble after each analysis step. It is a simple
     * way of increasing the variances/covariances with the effect of making KF corrections stronger.
     * The boolean variable turns it on or off, and the inflation factor controls how much deviations from the
     * ensemble mean are multiplied at each analysis step.
    */
    public boolean ensembleInflation = false;
    public double ensembleInflationFactor = 1.05;

    public double scaleAllPerturb = 1.0;

    // EOF based perturbations:
    public boolean eofPerturbations = true;
    public String eofPerturbationFile = "eof_perturb.nc";
    public double eofSurfScaleFactor = scaleAllPerturb*0.05;
    public double eofReductionRate = 0.03;

    public double ambientO2Std = 0.*scaleAllPerturb*0.25, ambientO2Beta = 0.2*0.05;

    public double currentStd = scaleAllPerturb*0.03, currentBeta = 0.2*0.05;
    public double o2ConsStd = scaleAllPerturb*0.15, o2ConsBeta = 0.2*0.05;

    // Parameter estimation:
    public int nPar = 1;
    double parStd0 = 0.002;
    double par4Std = 0.25;
    public double[] parStd = new double[] {parStd0};//, parStd0, parStd0};

    public double[] parDepths = new double[] {5, 10, 15};


    public boolean isDropOutActive(double timeS) {
        return correctionsDropOut && (timeS>correctionsDropOutStartS) && (timeS< correctionsDropOutEndS);
    }
}
