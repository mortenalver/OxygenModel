package fishmodel.hydraulics;

import ucar.ma2.ArrayDouble;
import ucar.ma2.ArrayFloat;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.NetcdfFile;
import ucar.nc2.Variable;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;


public class StationaryTankFlow {

    double[][][][] hydro;


    private float[] xc, yc, zc;
    double dxyFluent, dzFluent;
    double cOffsetX, cOffsetY;

    public StationaryTankFlow(String filepath, int[] cageDims, double dxy, double dz) {
        // First read the entire 3D field from the NetCDF file:
        float[][][] cfieldsU=null, cfieldsV=null, cfieldsW=null;
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
        dxyFluent = xc[1]-xc[0];
        dzFluent = zc[1]-zc[0];
        double centerX = 0.5*(xc[0]+ xc[xc.length-1]);
        double centerY = 0.5*(yc[0]+ yc[yc.length-1]);
        //System.out.println("CenterX = "+centerX+", CenterY = "+centerY);
        // Domain dimensions:
        double[] domDims = new double[] {cageDims[0]*dxy, cageDims[1]*dxy, cageDims[2]*dz};
        double domCenterX = domDims[0]/2., domCenterY = domDims[1]/2.;
        cOffsetX = domCenterX-centerX;
        cOffsetY = domCenterY-centerY;
        //System.out.println("Center offset X = "+cOffsetX);
        //System.out.println("Center offset Y = "+cOffsetY);
        // Define flow field variable:
        hydro = new double[cageDims[0]+1][cageDims[1]+1][cageDims[2]+1][3];


        for (int i = 0; i < cageDims[0] + 1; i++) {
            double xvalU = i * dxy - cOffsetX;
            double xvalV = (i + 0.5) * dxy - cOffsetX;
            //int xcoordU = (int) Math.round(xvalU / dxyFluent);
            //int xcoordV = (int) Math.round(xvalV / dxyFluent);
            //System.out.println("xValU="+xvalU+"; xccordU="+xcoordU);
            //System.out.println("xValV="+xvalV+"; xccordV="+xcoordV);
            for (int j = 0; j < cageDims[1] + 1; j++) {
                double yvalU = (j + 0.5) * dxy - cOffsetY;
                double yvalV = j * dxy - cOffsetY;
                //int ycoordU = (int) Math.round(yvalU / dxyFluent);
                //int ycoordV = (int) Math.round(yvalV / dxyFluent);
                //System.out.println("yValU="+yvalU+"; yccordU="+ycoordU);
                //System.out.println("yValV="+yvalV+"; yccordV="+ycoordV);
                for (int k = 0; k < cageDims[2] + 1; k++) {
                    double zvalU = zc[zc.length - 1] - (k + 0.5) * dxy;
                    double zvalW = zc[zc.length - 1] - k * dxy;
                    //System.out.println("z = "+zc[zc.length - 1]+" - "+k+" * "+dxy);
                    //int zcoordU = (int) Math.round(zvalU / dzFluent);
                    //int zcoordW = (int) Math.round(zvalW / dzFluent);
                    //System.out.println("zValU="+zvalU+"; zccordU="+zcoordU);
                    //System.out.println("zValW="+zvalW+"; zccordW="+zcoordW);

                    // U,V,W speeds:
                    hydro[i][j][k][0] = getValueU(cfieldsU, xvalU, yvalU, zvalU, dxy, dz);
                    hydro[i][j][k][1] = getValueV(cfieldsV, xvalV, yvalV, zvalU, dxy, dz);
                    boolean message = false;
                    //if (i==44 && j==19 && k==9) message = true;
                    if (k>0)
                        hydro[i][j][k][2] = -getValueW(cfieldsW, xvalV, yvalU, zvalW, dxy, dz, message);

                    /*if (Math.abs(hydro[i][j][k][2]) > 0.5) {
                        System.out.println("High val at: ("+i+", "+j+", "+k+"): "+hydro[i][j][k][2]);
                        System.out.println("Position: "+xvalV+", "+yvalU+", "+zvalW);
                    }*/
                }
            }
        }

        try {
            double imbalance = 0;

            FileWriter fw = new FileWriter(new File("balances.csv"));

            for (int k = 0; k < cageDims[2]; k++) {
                FileWriter fw2 = new FileWriter(new File("w_plane_"+k+".csv"));
                for (int i = 0; i < cageDims[0]; i++) {
                    for (int j = 0; j < cageDims[1]; j++) {

                        fw2.write(""+hydro[i][j][k][2]);
                        if (j<cageDims[1]-1)
                            fw2.write("\t");

                        double balance = hydro[i][j][k][0] - hydro[i + 1][j][k][0]
                                + hydro[i][j][k][1] - hydro[i][j + 1][k][1]
                                + hydro[i][j][k][2] - hydro[i][j][k + 1][2];
                        //System.out.println("Balance: "+balance);
                        imbalance += balance;

                        /*if (i==36 && j==36 && k==6) {
                            System.out.println("i="+i+", j="+j+", k="+k+", Balance: "+balance);
                            System.out.println("   "+hydro[i][j][k][0]+"   "+hydro[i + 1][j][k][0]+"   "
                                    + hydro[i][j][k][1]+"   "+hydro[i][j + 1][k][1]+"   "
                                    + hydro[i][j][k][2]+"   "+hydro[i][j][k + 1][2]);
                        }*/

                        fw.write(balance+"\r\n");

                    }
                    fw2.write("\r\n");
                }
                fw2.close();
            }
            fw.close();

            System.out.println("Total imbalance: "+imbalance);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

    }

    public int[] getCellForCoordinates(double[] position, double dxy, double dz) {
        int[] inds = new int[3];
        inds[0] = (int)Math.floor((position[0] + cOffsetX)/dxy);
        inds[1] = (int)Math.floor((position[1] + cOffsetY)/dxy);
        inds[2] = (int) Math.floor((zc[zc.length-1] - position[2])/dz);
        return inds;
    }

    public double[][][][] getFlowField() {
        return hydro;
    }

    private float getValueU(float[][][] a, double x, double y, double z, double dxy, double dz) {
        // All points should come from a single x coordinate:
        int xi = (int) Math.round(x / dxyFluent);
        // y and z coordinates should come from a range of +/- the O2 model grid size:
        int yi_min = (int) Math.round((y-0.5*dxy) / dxyFluent);
        int yi_max = (int) Math.round((y+0.5*dxy) / dxyFluent);
        int zi_min = (int) Math.round((z-0.5*dz) / dzFluent);
        int zi_max = (int) Math.round((z+0.5*dz) / dzFluent);

        int nval = 0;
        float sum = 0f;
        for (int i=yi_min; i<=yi_max; i++) {
            for (int j=zi_min; j<=zi_max; j++) {
                float vv = getValueN(a, xi, i, j);
                if (!Float.isNaN(vv)) {
                    sum += vv;
                    nval++;
                }
            }
        }
        if (nval > 0)
            return sum/((float)nval);
        else return 0f;
    }

    private float getValueV(float[][][] a, double x, double y, double z, double dxy, double dz) {
        // All points should come from a single x coordinate:
        int yi = (int) Math.round(y / dxyFluent);
        // y and z coordinates should come from a range of +/- the O2 model grid size:
        int xi_min = (int) Math.round((x-0.5*dxy) / dxyFluent);
        int xi_max = (int) Math.round((x+0.5*dxy) / dxyFluent);
        int zi_min = (int) Math.round((z-0.5*dz) / dzFluent);
        int zi_max = (int) Math.round((z+0.5*dz) / dzFluent);
        //System.out.println("x="+x+", y="+y+", z="+z);
        //System.out.print("xi: "+xi_min+" - "+xi_max+"\t");
        //System.out.print("yi: "+yi+"\t");
        //System.out.println("zi: "+zi_min+" - "+zi_max);

        int nval = 0;
        float sum = 0f;
        for (int i=xi_min; i<=xi_max; i++) {
            for (int j=zi_min; j<=zi_max; j++) {
                float vv = getValueN(a, i, yi, j);
                if (!Float.isNaN(vv)) {
                    sum += vv;
                    nval++;
                }
            }
        }
        if (nval > 0)
            return sum/((float)nval);
        else return 0f;
    }

    private float getValueW(float[][][] a, double x, double y, double z, double dxy, double dz, boolean message) {
        // All points should come from a single x coordinate:
        int zi = (int) Math.round(z / dzFluent);
        // y and z coordinates should come from a range of +/- the O2 model grid size:
        int xi_min = (int) Math.round((x-0.5*dxy) / dxyFluent);
        int xi_max = (int) Math.round((x+0.5*dxy) / dxyFluent);
        int yi_min = (int) Math.round((y-0.5*dxy) / dxyFluent);
        int yi_max = (int) Math.round((y+0.5*dxy) / dxyFluent);
        //System.out.println("x="+x+", y="+y+", z="+z);
        //System.out.print("xi: "+xi_min+" - "+xi_max+"\t");
        //System.out.println("yi: "+yi_min+" - "+yi_max);
        //System.out.print("zi: "+zi+"\t");

        int nval = 0;
        float sum = 0f;
        for (int i=xi_min; i<=xi_max; i++) {
            for (int j=yi_min; j<=yi_max; j++) {
                float vv = getValueN(a, i, j, zi);

                if (message)
                    System.out.println("vv = "+vv+", x="+x+", y="+y+", zi="+zi);

                if (!Float.isNaN(vv)) {
                    sum += vv;
                    nval++;
                }
            }
        }
        if (nval > 0)
            return sum/((float)nval);
        else return 0f;
    }

    private float getValueN(float[][][] a, int i, int j, int k) {
        if ((i<0) || (i>=a.length) || (j<0) || (j>=a[0].length) || (k<0) || (k>=a[0][0].length))
            return Float.NaN;
        return a[i][j][k];

    }

    private float getValue(float[][][] a, int i, int j, int k) {
        if ((i<0) || (i>=a.length) || (j<0) || (j>=a[0].length) || (k<0) || (k>=a[0][0].length))
            return 0.f;
        float val = a[i][j][k];
        if (Float.isNaN(val))
            return 0f;
        else return val;
    }


}
