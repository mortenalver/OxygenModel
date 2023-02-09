package fishmodel.hydraulics;

import ucar.ma2.ArrayFloat;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.NetcdfFile;
import ucar.nc2.Variable;

import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: malv
 * Date: 19.07.12
 * Time: 08:31
 * To change this template use File | Settings | File Templates.
 */
public class SimpleTankHydraulics {

    public static double[][][][] getUniformHydraulicField(int[] dims, double[] components) {
        double[][][][] field = new double[dims[0]+1][dims[1]+1][dims[2]+1][3];
        for (int k=0; k<field[0][0].length; k++)
            for (int i=0; i<field.length; i++)
                for (int j=0; j<field[0].length; j++) {
                    field[i][j][k][0] = components[0];
                    field[i][j][k][1] = components[1];
                    field[i][j][k][2] = components[2];
                }
        return field;
    }

    public static double[][][][] getCircularHydraulicField(double t, int[] dims, double maxspeed) {
        double[][][][] field = new double[dims[0]+1][dims[1]+1][dims[2]+1][3];
        double[] center = new double[3];
        for (int i=0; i<3; i++)
            center[i] = ((double)(dims[i]-1))/2.0;
        double radius = ((double)dims[0])/2.0;
        for (int k=0; k<field[0][0].length; k++)
            for (int i=0; i<field.length; i++)
                for (int j=0; j<field[0].length; j++) {
                    double distX = (double)i - 0.5 - center[0],
                            distY = (double)j - 0.5 - center[1];
                    double rpos = Math.sqrt(distX*distX + distY*distY);
                    if (rpos == 0) {
                        field[i][j][k][0] = 0;
                        field[i][j][k][1] = 0;
                        field[i][j][k][2] = 0;
                    }
                    else {
                        double speed = Math.max(0, maxspeed*(Math.cos(0.5*Math.PI*((rpos - 0.5*radius)*2/radius))));
                        //field[i][j][k][0] = speed;
                        //field[i][j][k][1] = rpos;
                        field[i][j][k][0] = -(speed/rpos)*distY;
                        field[i][j][k][1] = (speed/rpos)*distX;
                        field[i][j][k][2] = 0;
                    }
                }

        return field;
    }

    public static double[][][][] getNetCDFHydraulicField(String filepath, int depthLayers, double multiplier) {
        try {
            NetcdfFile ncfile = NetcdfFile.open(filepath);
            Variable u = ncfile.findVariable("u"), v = ncfile.findVariable("v");
            int[] shape = u.getShape();
            ArrayFloat.D2 udata = (ArrayFloat.D2)u.read(new int[] {0,0}, shape);
            ArrayFloat.D2 vdata = (ArrayFloat.D2)v.read(new int[] {0,0}, shape);
            double[][][][] hydro = new double[shape[0]][shape[1]][depthLayers+1][3];
            for (int i=0; i<shape[0]; i++)
                for (int j=0; j<shape[1]; j++) {
                    double uv = udata.get(i,j);
                    double vv = vdata.get(i,j);
                    for (int k=0; k<depthLayers; k++) {
                        hydro[i][j][k][0] = multiplier*uv;
                        hydro[i][j][k][1] = multiplier*vv;
                        hydro[i][j][k][2] = 0;
                    }
                }
            return hydro;
        } catch (IOException ex) {
            ex.printStackTrace();
            return null;
        } catch (InvalidRangeException e) {
            e.printStackTrace();
            return null;
        }
    }

    public static double[][][][] getProfileHydraulicField(int[] dims, double[][] profile) {
        double[][][][] field = new double[dims[0]+1][dims[1]+1][dims[2]+1][3];
        for (int k=0; k<field[0][0].length; k++)
            for (int i=0; i<field.length; i++)
                for (int j=0; j<field[0].length; j++) {
                    field[i][j][k][0] = profile[k][0];
                    field[i][j][k][1] = profile[k][1];
                    field[i][j][k][2] = profile[k][2];
                }
        return field;
    }

    public static void getProfileHydraulicField(double[][][][] field, int[] dims, double[][] profile) {
        for (int k=0; k<field[0][0].length; k++)
            for (int i=0; i<field.length; i++)
                for (int j=0; j<field[0].length; j++) {
                    field[i][j][k][0] = profile[k][0];
                    field[i][j][k][1] = profile[k][1];
                    field[i][j][k][2] = profile[k][2];
                }
    }
}
