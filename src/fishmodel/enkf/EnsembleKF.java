package fishmodel.enkf;

import fishmodel.Measurements;
import fishmodel.sim.InputDataNetcdf;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import save.SaveForEnKF;
import ucar.nc2.NetcdfFileWriteable;

import java.io.IOException;
import java.util.Random;

public class EnsembleKF {

    Random rand = new Random();
    boolean first = true;
    String filePrefix;
    NetcdfFileWriteable enKFFile = null;
    private int[] cageDims;
    private double dxy;
    Measurements.MeasurementSet measSet = null;
    DMatrixRMaj M = null;
    DMatrixRMaj M_allsensors = null;
    int callCount = 0;
    final int saveInterval = 10; // Number of calls between each time we save data to file

    EnOI enOI = null; // Reference to EnOI object. Instantiated if we are using a hybrid algorithm
    DMatrixRMaj K_enOI = null; // Reference to EnOI Kalman gain. Instantiated if we are using a hybrid algorithm

    public EnsembleKF(String filePrefix, int[] cageDims, int nPar, double dxy, Measurements.MeasurementSet measSet) {
        this.filePrefix = filePrefix;
        this.cageDims = cageDims;
        this.dxy = dxy;
        this.measSet = measSet;
        M = Measurements.getMeasurementModelBjoroya(cageDims, nPar, measSet, false);
        M_allsensors = Measurements.getMeasurementModelBjoroya(cageDims, nPar, measSet, true);

    }

    public double[][] doAnalysis(double t, double[][] dX, AssimSettings as,
                                 InputDataNetcdf inData) {

        int N = (as.useTwin ? dX.length-1: dX.length); // Set correct N corresponding to the ensemble X.
        //System.out.println("N = "+N);
        String file = filePrefix+"_ens.nc";
        boolean savingThisStep = false;
        if (first) {
            first = false;
            enKFFile = SaveForEnKF.initializeEnKFFile(file, dX[0].length, N, M.getNumRows(), M_allsensors.getNumRows(),
                    as.nPar, "seconds", cageDims, dxy, measSet);
            savingThisStep = true;

            if (as.hybrid_EnKF_ENOI) {
                enOI = new EnOI(filePrefix+"_enoi", as, cageDims, dxy, measSet);
                enOI.calculateKalmanGain(t, as);
                K_enOI = enOI.getKalmanGain();
                System.out.println("K_enOI: "+K_enOI.getNumRows()+" x "+K_enOI.getNumCols());

            }
        } else {
            if (++callCount == saveInterval) {
                enKFFile = SaveForEnKF.openFile(file);
                savingThisStep = true;
                callCount = 0;
            }
        }

        // Ensemble state:
        DMatrixRMaj X = CommonOps_DDRM.transpose(new DMatrixRMaj(dX), null);
        DMatrixRMaj X_twin = null;

        // Measurements: There are two main modes of operation. If "includeTwin" is True, we are obtaining measuremens
        // from a "twin" model running parallel to the ensemble. The set of states to measure is decided
        // by the measurement model provided by the measurements module.
        // If "useTwin" is False, we are using measurements loaded from file by the simInputsNetcdf module.
        DMatrixRMaj D=null, D_exact=null, D_all_exact = null;

        if (as.useTwin) {
            // Extract the first N-1 members as the ensemble, and the last as the twin:
            DMatrixRMaj X_ens = new DMatrixRMaj(X.getNumRows(),X.getNumCols()-1);
            CommonOps_DDRM.extract(X, 0, X.getNumRows(), 0,X.getNumCols()-1, X_ens);
            X_twin = new DMatrixRMaj(X.getNumRows(),1);
            CommonOps_DDRM.extract(X, 0, X.getNumRows(), X.getNumCols()-1,X.getNumCols(), X_twin);
            X = X_ens; // X now consists of the non-twin members
            // Get measurements from the twin:
            DMatrixRMaj d = CommonOps_DDRM.mult(M, X_twin, null);
            // Set up matrix with perturbed measurements for all ensemble members:
            DMatrixRMaj e_n = new DMatrixRMaj(1, X.getNumCols());
            e_n.fill(1.0);
            D = CommonOps_DDRM.mult(d, e_n, null);
            D_exact = D.copy(); // Store the unperturbed measurement matrix
            for (int i=0; i<D.getNumRows(); i++)
                for (int j=0; j<D.getNumCols(); j++)
                    D.add(i,j, measSet.std*rand.nextGaussian()); // Add gaussian noise according to measurement uncertainty
            // Get a record of all measurements (not just assimilated sensors):
            d = CommonOps_DDRM.mult(M_allsensors, X_twin, null);
            D_all_exact = CommonOps_DDRM.mult(d, e_n, null);


        }
        else {
            // Keep the entire X matrix, and let X_twin stay null. Get measurements from the input data object:
            double[] meas = inData.getO2Measurements();
            // Get a list of which sensors to use for assimilation:
            int[] activeSensors = Measurements.getSensorsToAssimilateBjoroya();

            DMatrixRMaj d = new DMatrixRMaj(activeSensors.length, 1);
            for (int i = 0; i < activeSensors.length; i++) {
                d.set(i, 0, meas[activeSensors[i]]);
            }
            // Set up matrix with perturbed measurements for all ensemble members:
            DMatrixRMaj e_n = new DMatrixRMaj(1, X.getNumCols());
            e_n.fill(1.0);
            D = CommonOps_DDRM.mult(d, e_n, null);
            D_exact = D.copy(); // Store the unperturbed measurement matrix
            for (int i=0; i<D.getNumRows(); i++)
                for (int j=0; j<D.getNumCols(); j++)
                    D.add(i,j, measSet.std*rand.nextGaussian()); // Add gaussian noise according to measurement uncertainty
            //System.out.println("D:"); System.out.println(D);

            // Get a record of all measurements (not just assimilated sensors):
            d = new DMatrixRMaj(meas.length, 1);
            for (int i = 0; i < meas.length; i++) {
                d.set(i, 0, meas[i]);
            }
            D_all_exact = CommonOps_DDRM.mult(d, e_n, null);
        }


        double N_d = (double)N;
        double N_1 = Math.max((double)N-1., 1.);
        int numel = cageDims[0]*cageDims[1]*cageDims[2]; // The number of cells in the field. May be lower than the
                // number of elements in the state vector, if parameters are being estimated.

        // Set up R matrix:
        double[] rval = new double[M.getNumRows()];
        for (int i = 0; i < rval.length; i++) {
            rval[i] = measSet.std*measSet.std;
        }
        DMatrixRMaj R = CommonOps_DDRM.diag(rval);

        // Compute mean state:
        DMatrixRMaj X_mean = CommonOps_DDRM.sumRows(X, null);
        CommonOps_DDRM.scale(1./N_d, X_mean);
        //System.out.println("X_mean: "+X_mean.getNumRows()+"x"+X_mean.getNumCols());
        DMatrixRMaj e_n = new DMatrixRMaj(1, X.getNumCols());
        e_n.fill(1.0);
        DMatrixRMaj E_X = CommonOps_DDRM.mult(X_mean, e_n, null);

        // Compute deviations from mean state:
        DMatrixRMaj theta = CommonOps_DDRM.subtract(X, E_X, null);
        //System.out.println("theta: "+theta.getNumRows()+"x"+theta.getNumCols());

        // Compute intermediary matrixes:
        DMatrixRMaj MX = CommonOps_DDRM.mult(M, X, null);
        DMatrixRMaj MX_allsensors = CommonOps_DDRM.mult(M_allsensors, X, null);
        //System.out.println("M: "+M.getNumRows()+"x"+M.getNumCols());
        //System.out.println("MX: "+MX.getNumRows()+"x"+MX.getNumCols());
        DMatrixRMaj omega = CommonOps_DDRM.mult(M, theta, null);
        //System.out.println("omega: "+omega.getNumRows()+"x"+omega.getNumCols());

        //phi = (1./N1)*(omega @ omega.T) + R
        DMatrixRMaj p1 = CommonOps_DDRM.mult(omega, CommonOps_DDRM.transpose(omega, null), null);
        CommonOps_DDRM.scale(1./N_1, p1);
        DMatrixRMaj phi_inv = CommonOps_DDRM.add(p1, R, null);
        CommonOps_DDRM.invert(phi_inv);
        //System.out.println("phi_inv: "+phi_inv.getNumRows()+"x"+phi_inv.getNumCols());

        System.out.println("D: "+D.getNumRows()+"x"+D.getNumCols());
        DMatrixRMaj deviations = CommonOps_DDRM.subtract(D, MX, null);
        DMatrixRMaj dev_exact = CommonOps_DDRM.subtract(D_exact, MX, null);
        DMatrixRMaj dev_all_exact = CommonOps_DDRM.subtract(D_all_exact, MX_allsensors, null);
        //System.out.println("Deviations:");
        //System.out.println(deviations);
        //System.out.println("Dev exact:");
        //System.out.println(dev_exact);

        DMatrixRMaj Kloc1 = getLocalizationMatrix(cageDims, M, numel, as.locDist, as.locZMultiplier);

        // Calculate Kalman gain:
        DMatrixRMaj K = CommonOps_DDRM.elementMult(Kloc1, CommonOps_DDRM.mult(CommonOps_DDRM.mult(theta, CommonOps_DDRM.transpose(omega,
                null), null), phi_inv, null), null);
        System.out.println("K: "+K.getNumRows()+" x "+K.getNumCols());
        CommonOps_DDRM.scale(1./N_1, K);

        DMatrixRMaj corrections = null;
        // Check if we are running a hybrid setup:
        if (as.hybrid_EnKF_ENOI) {
            DMatrixRMaj K_EnOI_now = null;
            // Calculate a weighted average of the EnKF and EnOI kalman gains.
            if (as.nPar > 0) {
                // We need to copy K elements for parameter adjustments from the EnKF K since EnOI doesn't include parameters:
                DMatrixRMaj parRows = new DMatrixRMaj(as.nPar, K.getNumCols());
                CommonOps_DDRM.extract(K, K.getNumRows()-as.nPar, K.getNumRows(), 0, K.getNumCols(), parRows);
                System.out.println("parRows: "+parRows.getNumRows()+" x "+parRows.getNumCols());
                K_EnOI_now = CommonOps_DDRM.concatRowsMulti(K_enOI, parRows);
                System.out.println("K_EnOI_now: "+K_EnOI_now.getNumRows()+" x "+K_EnOI_now.getNumCols());
            } else {
                K_EnOI_now = K_enOI;
            }
            CommonOps_DDRM.scale(as.hybrid_ENOI_weight, K_EnOI_now);
            DMatrixRMaj K_copy = new DMatrixRMaj(K);
            CommonOps_DDRM.scale(1.0-as.hybrid_ENOI_weight, K_copy);
            DMatrixRMaj K_weighted = CommonOps_DDRM.add(K_copy, K_EnOI_now, null);
            System.out.println("K_weighted: "+K_weighted.getNumRows()+" x "+K_weighted.getNumCols());
            corrections = CommonOps_DDRM.mult(K_weighted, deviations, null);
        }
        else {
            // Not hybrid, so just go ahead with the K we computed:
            corrections = CommonOps_DDRM.mult(K, deviations, null);
        }
        //System.out.println("Corrections: "+corrections.getNumRows()+"x"+corrections.getNumCols());

        DMatrixRMaj X_a = CommonOps_DDRM.add(X, corrections, null);

        // Ensemble inflation:
        if (as.ensembleInflation) {
            // Compute mean of analysis:
            CommonOps_DDRM.sumRows(X_a, X_mean);
            CommonOps_DDRM.scale(1./N_d, X_mean);
            // Expand mean into N identical columns:
            CommonOps_DDRM.mult(X_mean, e_n, E_X);
            // Compute X_a - E_X:
            DMatrixRMaj X_a_minus_E_X = CommonOps_DDRM.subtract(X_a, E_X, null);
            // Scale it by the inflation factor:
            CommonOps_DDRM.scale(as.ensembleInflationFactor, X_a_minus_E_X);
            // Compute the inflated ensemble state and put it into X_a:
            CommonOps_DDRM.add(X_a_minus_E_X, E_X, X_a);
        }

        // Save ensemble state:
        if (savingThisStep) {
            SaveForEnKF.saveEnKFVariables(enKFFile, t, m2A(X), m2A(X_a),
                    (X_twin != null ? m2A_1d(X_twin) : null), m2A(M),
                    m2A(deviations), m2A(dev_exact), m2A(dev_all_exact), m2A(K), m2A(Kloc1));
            try {
                enKFFile.close();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
        return m2A_transpose(X_a);
    }

    public DMatrixRMaj getLocalizationMatrix(int[] dims, DMatrixRMaj M, int numel, double locDist, double locZMultiplier) {
        DMatrixRMaj Xloc1 = new DMatrixRMaj(M.getNumCols(), M.getNumRows());

        for (int i=0; i<M.getNumRows(); i++) { // Loop over measurements
            // Find index of this measurement:
            int mInd = 0;
            while (M.get(i, mInd) == 0)
                mInd = mInd + 1;
            double[] mCoord = Util.getStateCoords(mInd, dims);

            // Set values for Xloc1:
            for (int j = 0; j < M.getNumCols(); j++) { // Loop over state variables
                if (j < numel) {
                    double[] sCoord = Util.getStateCoords(j, dims);
                    double distance = Math.sqrt(Math.pow(mCoord[0] - sCoord[0], 2) + Math.pow(mCoord[1] - sCoord[1], 2) +
                            Math.pow(locZMultiplier * (mCoord[2] - sCoord[2]), 2));
                    Xloc1.set(j, i, localizationValue(distance, locDist));
                } else {
                    // No localization for parameter values:
                    Xloc1.set(j, i, 1.);
                }
            }
        }
        return Xloc1;
    }

    private double localizationValue(double dist, double locDist) {
        // Localization function used in Daniel's EnKF:
        double lval;
        double c = Math.sqrt(10./3.)*locDist;
        double frac = dist/c;
        if (dist<c) {
            lval = -0.25 * Math.pow(frac, 5.) + 0.5 * Math.pow(frac, 4.) + (5. / 8.) * Math.pow(frac, 3.) - (5. / 3.) * Math.pow(frac, 2.) + 1.;
        }
        else if (dist <2 * c) {
            lval = (1. / 12.) * Math.pow(frac, 5.) - 0.5 * Math.pow(frac, 4.) + (5. / 8.) * Math.pow(frac, 3.) + (5. / 3.) * Math.pow(frac, 2.) -
                5. * frac + 4. - (2. / 3.) / frac;
        }
        else
            lval = 0.;

        return lval;
    }

    public double[][] m2A(DMatrixRMaj matrix) {
        double[][] array = new double[matrix.getNumRows()][matrix.getNumCols()];
        for (int r = 0; r < array.length; r++) {
            for (int c = 0; c < array[r].length; c++) {
                array[r][c] = matrix.get(r, c);
            }
        }
        return array;
    }

    public double[][] m2A_transpose(DMatrixRMaj matrix) {
        double[][] array = new double[matrix.getNumCols()][matrix.getNumRows()];
        for (int r = 0; r < array.length; r++) {
            for (int c = 0; c < array[r].length; c++) {
                array[r][c] = matrix.get(c, r);
            }
        }
        return array;
    }

    public double[] m2A_1d(DMatrixRMaj matrix) {
        double[] array = new double[matrix.getNumRows()*matrix.getNumCols()];
        int index = 0;
        for (int r = 0; r < matrix.getNumRows(); r++) {
            for (int c = 0; c < matrix.getNumCols(); c++) {
                array[index++] = matrix.get(r, c);
            }
        }
        return array;
    }

    public static void main(String[] args) {
        DMatrixRMaj X = new DMatrixRMaj(2,3);
        X.set(0,0,1); X.set(0,1,6); X.set(0,2,3);
        X.set(1,0,2); X.set(1,1,3); X.set(1,2,4);

        DMatrixRMaj X2 = new DMatrixRMaj(X.getNumRows(),X.getNumCols()-1);
        CommonOps_DDRM.extract(X, 0, X.getNumRows(), 0,X.getNumCols()-1, X2);
        DMatrixRMaj X3 = new DMatrixRMaj(X.getNumRows(),1);
        CommonOps_DDRM.extract(X, 0, X.getNumRows(), X.getNumCols()-1,X.getNumCols(), X3);

        System.out.println(X);
        System.out.println(X2);
        System.out.println(X3);



    }
}
