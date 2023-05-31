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
        String file = filePrefix+"_ens.nc";
        if (first) {
            first = false;
            enKFFile = SaveForEnKF.initializeEnKFFile(file, dX[0].length, dX.length, 1, "seconds");

        } else {
            enKFFile = SaveForEnKF.openFile(file);
        }
        SaveForEnKF.saveEnKFVariables(enKFFile, t, dX);

        // Ensemble state:
        DMatrixRMaj X = CommonOps_DDRM.transpose(new DMatrixRMaj(dX), null);

        // Measurements: There are two main modes of operation. If "includeTwin" is True, we are obtaining measuremens
        // from a "twin" model running parallel to the ensemble. The set of states to measure is decided
        // by the measurement model provided by the measurements module.
        // If "useTwin" is False, we are using measurements loaded from file by the simInputsNetcdf module.
        if (useTwin) {
            // Extract the first N-1 members as the ensemble, and the last as the twin:
            DMatrixRMaj X_ens = new DMatrixRMaj(X.getNumRows(),X.getNumCols()-1);
            CommonOps_DDRM.extract(X, 0, X.getNumRows(), 0,X.getNumCols()-1, X_ens);
            DMatrixRMaj X_twin = new DMatrixRMaj(X.getNumRows(),1);
            CommonOps_DDRM.extract(X, 0, X.getNumRows(), X.getNumCols()-1,X.getNumCols(), X_twin);
            X = X_ens; // X now consists of the non-twin members
            // Get measurements from the twin:
            DMatrixRMaj d = CommonOps_DDRM.mult(M, X_twin, null);
            System.out.println("Measurement:"); System.out.println(d);
            // Set up matrix with perturbed measurements for all ensemble members:
            DMatrixRMaj e_n = new DMatrixRMaj(1, X.getNumCols());
            e_n.fill(1.0);
            DMatrixRMaj D = CommonOps_DDRM.mult(d, e_n, null);
            for (int i=0; i<D.getNumRows(); i++)
                for (int j=0; j<D.getNumCols(); j++)
                    D.add(i,j, measSet.std*rand.nextGaussian());
            System.out.println("D:"); System.out.println(D);

        } {
            // TODO: handle non-twin case
        }

        int N = X.getNumCols(); // Set correct N corresponding to the ensemble X.

        //CommonOps_DDRM.extract();
        //DMatrixRMaj mXX = CommonOps_DDRM.mult(mX, CommonOps_DDRM.transpose(mX, null), null);
        //System.out.println("XX: "+mXX.getNumRows()+" x "+mXX.getNumCols());
        //System.out.println(mXX);

        try {
            enKFFile.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
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
