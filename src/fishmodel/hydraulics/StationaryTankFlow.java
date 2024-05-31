package fishmodel.hydraulics;

import ucar.ma2.ArrayDouble;
import ucar.ma2.ArrayFloat;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.NetcdfFile;
import ucar.nc2.Variable;

import java.io.IOException;


public class StationaryTankFlow {

    public static void main(String[] args) {
        StationaryTankFlow cmf = new StationaryTankFlow("C:\\Users\\alver\\OneDrive - NTNU\\prosjekt\\AQUAEXCEL3\\particleTransport\\inrae_mod_max_3d_4.nc");

    }

    private float[][][] cfieldsU, cfieldsV, cfieldsW;
    private float[] xc, yc, zc;
    public StationaryTankFlow(String filepath) {
        // First read the entire 3D field from the NetCDF file:
        try {
            NetcdfFile ncfile = NetcdfFile.open(filepath);

            // Read grid coordinates:
            Variable xcVar = ncfile.findVariable("xc");
            Variable ycVar = ncfile.findVariable("yc");
            Variable zcVar = ncfile.findVariable("zc");
            int[] xcShape = xcVar.getShape();
            int[] ycShape = ycVar.getShape();
            int[] zcShape = zcVar.getShape();
            xc = new float[xcShape[0]];
            ArrayFloat.D1 valueXc = (ArrayFloat.D1)xcVar.read(new int[] {0}, new int[] {xcShape[0]});
            for (int i=0; i<xcShape[0]; i++) {
                xc[i] = valueXc.get(i);
            }
            yc = new float[ycShape[0]];
            ArrayFloat.D1 valueYc = (ArrayFloat.D1)ycVar.read(new int[] {0}, new int[] {ycShape[0]});
            for (int i=0; i<ycShape[0]; i++) {
                yc[i] = valueYc.get(i);
            }
            zc = new float[zcShape[0]];
            ArrayFloat.D1 valueZc = (ArrayFloat.D1)zcVar.read(new int[] {0}, new int[] {zcShape[0]});
            for (int i=0; i<zcShape[0]; i++) {
                zc[i] = valueZc.get(i);
            }

            Variable uVar = ncfile.findVariable("u");
            Variable vVar = ncfile.findVariable("v");
            Variable wVar = ncfile.findVariable("w");

            // Get shape of u variable:
            int[] shape = uVar.getShape();
            for (int i = 0; i < shape.length; i++) {
                int i1 = shape[i];
                System.out.println("Shape["+i+"]="+i1);
            }

            cfieldsU = new float[shape[0]][shape[1]][shape[2]];
            cfieldsV = new float[shape[0]][shape[1]][shape[2]];
            cfieldsW = new float[shape[0]][shape[1]][shape[2]];
            for (int k=0; k<shape[2]; k++) {
                ArrayFloat.D3 valueU = (ArrayFloat.D3)uVar.read(new int[] {0, 0, k}, new int[] {shape[0], shape[1], 1});
                ArrayFloat.D3 valueV = (ArrayFloat.D3)vVar.read(new int[] {0, 0, k}, new int[] {shape[0], shape[1], 1});
                ArrayFloat.D3 valueW = (ArrayFloat.D3)wVar.read(new int[] {0, 0, k}, new int[] {shape[0], shape[1], 1});
                for (int i = 0; i < shape[0]; i++) {
                    for (int j = 0; j < shape[1]; j++) {
                        cfieldsU[i][j][k] = valueU.get(i,j,0);
                        cfieldsV[i][j][k] = valueV.get(i,j,0);
                        cfieldsW[i][j][k] = valueW.get(i,j,0);
                    }
                }
                System.out.println("done reading k="+k);
            }

            ncfile.close();
        } catch (
                IOException ex) {
            ex.printStackTrace();

        } catch (InvalidRangeException e) {
            throw new RuntimeException(e);
        }

        // Then we need to extract numbers for our model grid. This requires us to line up the
        // model grids properly.
        float centerX = 0.5f*(xc[0]+ xc[xc.length-1]);
        float centerY = 0.5f*(yc[0]+ yc[yc.length-1]);
        System.out.println("CenterX = "+centerX+", CenterY = "+centerY);

    }





}
