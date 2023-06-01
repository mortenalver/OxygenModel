package fishmodel.enkf;

import fishmodel.Measurements;
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
    Measurements.MeasurementSet measSet = null;
    DMatrixRMaj M = null;
    public EnsembleKF(String filePrefix, int[] cageDims, Measurements.MeasurementSet measSet) {
        this.filePrefix = filePrefix;
        this.measSet = measSet;
        M = Measurements.getMeasurementModel(cageDims, measSet);
    }

    public void doAnalysis(double t, double[][] dX, boolean useTwin) {

        int N = (useTwin ? dX.length-1: dX.length); // Set correct N corresponding to the ensemble X.
        System.out.println("N = "+N);
        String file = filePrefix+"_ens.nc";
        if (first) {
            first = false;
            enKFFile = SaveForEnKF.initializeEnKFFile(file, dX[0].length, N, 1, "seconds");

        } else {
            enKFFile = SaveForEnKF.openFile(file);
        }

        // Ensemble state:
        DMatrixRMaj X = CommonOps_DDRM.transpose(new DMatrixRMaj(dX), null);
        DMatrixRMaj X_twin = null;

        // Measurements: There are two main modes of operation. If "includeTwin" is True, we are obtaining measuremens
        // from a "twin" model running parallel to the ensemble. The set of states to measure is decided
        // by the measurement model provided by the measurements module.
        // If "useTwin" is False, we are using measurements loaded from file by the simInputsNetcdf module.
        DMatrixRMaj D=null, D_exact=null;
        if (useTwin) {
            // Extract the first N-1 members as the ensemble, and the last as the twin:
            DMatrixRMaj X_ens = new DMatrixRMaj(X.getNumRows(),X.getNumCols()-1);
            CommonOps_DDRM.extract(X, 0, X.getNumRows(), 0,X.getNumCols()-1, X_ens);
            X_twin = new DMatrixRMaj(X.getNumRows(),1);
            CommonOps_DDRM.extract(X, 0, X.getNumRows(), X.getNumCols()-1,X.getNumCols(), X_twin);
            X = X_ens; // X now consists of the non-twin members
            // Get measurements from the twin:
            DMatrixRMaj d = CommonOps_DDRM.mult(M, X_twin, null);
            System.out.println("Measurement:"); System.out.println(d);
            // Set up matrix with perturbed measurements for all ensemble members:
            DMatrixRMaj e_n = new DMatrixRMaj(1, X.getNumCols());
            e_n.fill(1.0);
            D = CommonOps_DDRM.mult(d, e_n, null);
            D_exact = D.copy(); // Store the unperturbed measurement matrix
            for (int i=0; i<D.getNumRows(); i++)
                for (int j=0; j<D.getNumCols(); j++)
                    D.add(i,j, measSet.std*rand.nextGaussian()); // Add gaussian noise according to measurement uncertainty
            System.out.println("D:"); System.out.println(D);

        } {
            // TODO: handle non-twin case
        }


        double N_d = (double)N;
        double N_1 = Math.max((double)N-1., 1.);
        // TODO: Localization

        // Set up R matrix:
        double[] rval = new double[M.getNumRows()];
        for (int i = 0; i < rval.length; i++) {
            rval[i] = measSet.std*measSet.std;
        }
        DMatrixRMaj R = CommonOps_DDRM.diag(rval);

        // Compute mean state:
        DMatrixRMaj X_mean = CommonOps_DDRM.sumRows(X, null);
        CommonOps_DDRM.scale(1./N_d, X_mean);
        System.out.println("X_mean: "+X_mean.getNumRows()+"x"+X_mean.getNumCols());
        DMatrixRMaj e_n = new DMatrixRMaj(1, X.getNumCols());
        e_n.fill(1.0);
        DMatrixRMaj E_X = CommonOps_DDRM.mult(X_mean, e_n, null);

        // Compute deviations from mean state:
        DMatrixRMaj theta = CommonOps_DDRM.subtract(X, E_X, null);
        System.out.println("theta: "+theta.getNumRows()+"x"+theta.getNumCols());

        // Compute intermediary matrixes:
        DMatrixRMaj MX = CommonOps_DDRM.mult(M, X, null);
        System.out.println("MX: "+MX.getNumRows()+"x"+MX.getNumCols());
        DMatrixRMaj omega = CommonOps_DDRM.mult(M, theta, null);
        System.out.println("omega: "+omega.getNumRows()+"x"+omega.getNumCols());

        //phi = (1./N1)*(omega @ omega.T) + R
        DMatrixRMaj p1 = CommonOps_DDRM.mult(omega, CommonOps_DDRM.transpose(omega, null), null);
        CommonOps_DDRM.scale(1./N_1, p1);
        DMatrixRMaj phi_inv = CommonOps_DDRM.add(p1, R, null);
        CommonOps_DDRM.invert(phi_inv);
        System.out.println("phi_inv: "+phi_inv.getNumRows()+"x"+phi_inv.getNumCols());

        DMatrixRMaj deviations = CommonOps_DDRM.subtract(D, MX, null);
        DMatrixRMaj dev_exact = CommonOps_DDRM.subtract(D_exact, MX, null);
        //System.out.println("Deviations:");
        //System.out.println(deviations);
        //System.out.println("Dev exact:");
        //System.out.println(dev_exact);

        // TODO: localization

        // Calculate Kalman gain:
        DMatrixRMaj K = CommonOps_DDRM.mult(CommonOps_DDRM.mult(theta, CommonOps_DDRM.transpose(omega, null), null), phi_inv, null);
        CommonOps_DDRM.scale(1./N_1, K);

        DMatrixRMaj corrections = CommonOps_DDRM.mult(K, deviations, null);
        System.out.println("Corrections: "+corrections.getNumRows()+"x"+corrections.getNumCols());

        DMatrixRMaj X_a = CommonOps_DDRM.add(X, corrections, null);

        // TODO: ensemble inflation?


        // Save ensemble state:
        SaveForEnKF.saveEnKFVariables(enKFFile, t, m2A(X), m2A(X_a), m2A_1d(X_twin), m2A(M));


        try {
            enKFFile.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
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
