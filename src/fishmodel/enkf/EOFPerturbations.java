package fishmodel.enkf;

import ucar.ma2.ArrayDouble;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.Dimension;
import ucar.nc2.NetcdfFile;
import ucar.nc2.Variable;

import java.io.IOException;
import java.util.Random;

public class EOFPerturbations {

    final static String varname = "perturbations";
    private String filename;
    private double[] scaleFactors = null;
    private int nRep = -1;

    public EOFPerturbations(String filename) {
        this.filename = filename;

    }

    public double[] getScaleFactors(double dx, int nLayers, double surfVal, double reductionRate) {
        double[] res = new double[nLayers];
        for (int i=0; i<res.length; i++) {
            double depth = (i + 0.5)*dx;
            res[i] = surfVal*Math.exp(-depth*reductionRate);
        }
        return res;
    }

    public void setScaleFactors(double dx, int nLayers, double surfVal, double reductionRate) {
        this.scaleFactors = getScaleFactors(dx, nLayers, surfVal, reductionRate);
    }

    public void perturb3dField(double[][][] field) {
        double[][] rndField = getRandomPerturbationField();
        for (int i=0; i<field.length; i++)
            for (int j=0; j<field[i].length; j++)
                for (int k=0; k<field[i][j].length; k++)
                    field[i][j][k] += scaleFactors[k]*rndField[i][j];
    }

    /**
     * Load a random 2D field of perturbations from the given NetCDF file, scaled by the given scale factor
     * @return The random 2D field
     */
    public double[][] getRandomPerturbationField() {
        try {
            NetcdfFile ncfile = NetcdfFile.open(filename);
            if (nRep < 0) { // First time, need to figure out how many examples we can choose from:
                Dimension dim = ncfile.findDimension("rep");
                nRep = dim.getLength();
            }
            int repToLoad = (int)Math.floor(nRep*Math.random());
            Variable v = ncfile.findVariable(varname);
            int[] shape = v.getShape();
            /*for (int i = 0; i < shape.length; i++) {
                int i1 = shape[i];
                System.out.println(i1);
            }*/
            ArrayDouble.D3 vdata = (ArrayDouble.D3)v.read(new int[]{repToLoad, 0, 0},
                    new int[] {1, shape[1], shape[2]});
            ncfile.close();
            double[][] res = new double[shape[1]][shape[2]];
            for (int i=0; i<shape[1]; i++)
                for (int j=0; j<shape[2]; j++)
                    res[i][j] = vdata.get(0, i, j);
            return res;
        } catch (IOException e) {
            throw new RuntimeException(e);
        } catch (InvalidRangeException e) {
            throw new RuntimeException(e);
        }


    }

    public static void main(String[] args) {
        EOFPerturbations ep = new EOFPerturbations("data/eof_perturb.nc");

        double[] scaleF = ep.getScaleFactors(2, 10, 1, 0.03);
        for (int i = 0; i < scaleF.length; i++) {
            double v = scaleF[i];
            System.out.println(v);
        }

    }
}
