package fishmodel.enkf;

import fishmodel.Measurements;
import fishmodel.sim.InputDataNetcdf;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import save.SaveForEnOI;
import ucar.ma2.ArrayDouble;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.NetcdfFile;
import ucar.nc2.NetcdfFileWriteable;
import ucar.nc2.Variable;

import java.io.IOException;
import java.util.Random;

public class EnOI {

    Random rand = new Random();
    boolean first = true;
    String filePrefix;
    private int[] cageDims;
    private double dxy;
    Measurements.MeasurementSet measSet = null;

    NetcdfFileWriteable saveFile = null;
    DMatrixRMaj X = null;
    DMatrixRMaj M = null;
    DMatrixRMaj M_allsensors = null;
    DMatrixRMaj K = null;
    int callCount = 0;
    final int saveInterval = 10; // Number of calls between each time we save data to file
    public EnOI(String prefix, AssimSettings as, int[] cageDims, double dxy, Measurements.MeasurementSet measSet) {
        this.filePrefix = prefix;
        this.cageDims = cageDims;
        this.dxy = dxy;
        this.measSet = measSet;
        M = Measurements.getMeasurementModelBjoroya(cageDims, 0, measSet, false);
        M_allsensors = Measurements.getMeasurementModelBjoroya(cageDims, 0, measSet, true);

    }

    private DMatrixRMaj loadEnsembleFromFile(String path) {
        try {
            NetcdfFile ncfile = NetcdfFile.open(path);
            Variable o2 = ncfile.findVariable("o2");
            int[] shape = o2.getShape();
            System.out.println("Got variable and shape: "+shape[0]+" / "+shape[1]+" / "+shape[2]+" / "+shape[3]);
            ArrayDouble.D4 o2Data = (ArrayDouble.D4) o2.read(new int[]{0, 0, 0, 0}, shape);
            ncfile.close();

            // Create ensemble matrix:
            int nStates = cageDims[0]*cageDims[1]*cageDims[2];
            DMatrixRMaj myX = new DMatrixRMaj(nStates, shape[0]);

            // Go through time steps and gather state:
            for (int tstep=0; tstep<shape[0]; tstep++) {
                int index = 0;
                for (int i = 0; i < cageDims[0]; i++) {
                    for (int j = 0; j < cageDims[1]; j++) {
                        for (int k = 0; k < cageDims[2]; k++) {
                            myX.set(index, tstep, o2Data.get(tstep, k, j, i));
                            index++;
                        }
                    }
                }
            }

            return myX;
        } catch (IOException e) {
            throw new RuntimeException(e);
        } catch (InvalidRangeException e) {
            throw new RuntimeException(e);
        }
    }

    public void loadEnsemble(String[] path) {
        System.out.println("Loading EnOI ensemble...");
        DMatrixRMaj[] allX = new DMatrixRMaj[path.length];
        for (int i = 0; i < path.length; i++) {
            allX[i] = loadEnsembleFromFile(path[i]);
            //System.out.println("dim: "+allX[i].getNumRows()+" x "+allX[i].getNumCols());
            if (allX[i] == null) {
                System.out.println("Ensemble file could not be loaded: "+path[i]);
                System.exit(0);
            }
        }
        X = CommonOps_DDRM.concatColumnsMulti(allX);
        //System.out.println("Dim X: "+X.getNumRows()+" x "+X.getNumCols());
    }

    public double[][] doAnalysis(double t, double[] x_f_array, AssimSettings as,
                                 InputDataNetcdf inData) {

        String file = filePrefix+"_ens.nc";
        boolean savingThisStep = false;
        if (first) {
            first = false;
            savingThisStep = true;
            calculateKalmanGain(t, as);

        } else {
            if (++callCount == saveInterval) {
                saveFile = SaveForEnOI.openFile(file);
                savingThisStep = true;
                callCount = 0;
            }
        }

        // Instantiate forecast state vector:
        DMatrixRMaj x_f = new DMatrixRMaj(x_f_array);
        // Get measurements from the input data object:
        double[] meas = inData.getO2Measurements();
        // Get a list of which sensors to use for assimilation:
        int[] activeSensors = Measurements.getSensorsToAssimilateBjoroya();
        DMatrixRMaj d_all = new DMatrixRMaj(meas.length, 1);
        for (int i=0; i<meas.length; i++)
            d_all.set(i, 0, meas[i]);
        DMatrixRMaj d = new DMatrixRMaj(activeSensors.length, 1);
        for (int i = 0; i < activeSensors.length; i++) {
            d.set(i, 0, meas[activeSensors[i]]);
        }

        // Get model estimates based on the measurement model:
        DMatrixRMaj y_est = CommonOps_DDRM.mult(M, x_f, null);
        DMatrixRMaj y_est_all = CommonOps_DDRM.mult(M_allsensors, x_f, null);


        DMatrixRMaj dev = CommonOps_DDRM.subtract(d, y_est, null);
        DMatrixRMaj dev_all = CommonOps_DDRM.subtract(d_all, y_est_all, null);

        DMatrixRMaj x_a = CommonOps_DDRM.add(x_f, CommonOps_DDRM.mult(K, dev, null), null);


        if (savingThisStep) {

            SaveForEnOI.saveVariables(saveFile, t, m2A_1d(x_f), m2A_1d(x_a), m2A_1d(dev), m2A_1d(dev_all));
            try {
                saveFile.close();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
        return m2A_transpose(x_a);
    }

    public DMatrixRMaj getKalmanGain() {
        return K;
    }

    public void calculateKalmanGain(double t, AssimSettings as) {

        String file = filePrefix+"_ens.nc";

        int N=0,n=0, numel=cageDims[0]*cageDims[1]*cageDims[2];

        // Load ensemble and compute Kalman gain:
        loadEnsemble(as.enOIEnsembleFile);
        N = X.getNumCols();
        n = X.getNumRows();
        double N_d = (double)N;
        double N_1 = Math.max((double)N-1., 1.);


        // Localization matrix:
        DMatrixRMaj Kloc = getLocalizationMatrix(cageDims, M, numel, as.locDist, as.locZMultiplier);

        // Set up R matrix:
        double[] rval = new double[M.getNumRows()];
        for (int i = 0; i < rval.length; i++) {
            rval[i] = measSet.std*measSet.std;
        }
        DMatrixRMaj R = CommonOps_DDRM.diag(rval);

        // Compute mean ensemble state:
        DMatrixRMaj X_mean = CommonOps_DDRM.sumRows(X, null);
        CommonOps_DDRM.scale(1./N_d, X_mean);
        DMatrixRMaj e_n = new DMatrixRMaj(1, X.getNumCols());
        e_n.fill(1.0);
        DMatrixRMaj E_X = CommonOps_DDRM.mult(X_mean, e_n, null);
        // Compute ensemble anomalies; A in Oke et al. (2010):
        DMatrixRMaj theta = CommonOps_DDRM.subtract(X, E_X, null);
        // Multiply anomalies by EnOI alpha:
        CommonOps_DDRM.scale(as.enoiAlpha, theta);
        // Compute M multiplied by ensemble anomalies; HA in Oke:
        DMatrixRMaj omega = CommonOps_DDRM.mult(M, theta, null);

        DMatrixRMaj p1 = CommonOps_DDRM.mult(omega, CommonOps_DDRM.transpose(omega, null), null);
        DMatrixRMaj mR = new DMatrixRMaj(R);
        CommonOps_DDRM.scale(N_1, mR);
        DMatrixRMaj phi_inv = CommonOps_DDRM.add(p1, mR, null);
        CommonOps_DDRM.invert(phi_inv);


        // Calculate Kalman gain:
        K = CommonOps_DDRM.elementMult(Kloc, CommonOps_DDRM.mult(CommonOps_DDRM.mult(theta, CommonOps_DDRM.transpose(omega,
                null), null), phi_inv, null), null);

        saveFile = SaveForEnOI.initializeFile(file, n, N, M.getNumRows(), M_allsensors.getNumRows(),
                as.nPar, "seconds", cageDims, dxy, measSet);
        SaveForEnOI.saveStaticVariables(saveFile, t, m2A(X), m2A(Kloc), m2A(K));
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


}
