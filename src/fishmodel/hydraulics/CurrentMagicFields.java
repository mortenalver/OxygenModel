package fishmodel.hydraulics;

import ucar.ma2.ArrayDouble;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.NetcdfFile;
import ucar.nc2.Variable;

import java.io.IOException;


public class CurrentMagicFields {

    public static void main(String[] args) {
        CurrentMagicFields cmf = new CurrentMagicFields("C:/Users/alver/OneDrive - NTNU/prosjekt/O2_Bjørøya/currentmagic/currents_bjoroya.nc");

    }
    private double[] thetas;
    private double d_theta = 0;
    private double[][][] cfieldsU, cfieldsV;

    public CurrentMagicFields(String filepath) {
        try {
            NetcdfFile ncfile = NetcdfFile.open(filepath);
            Variable uVar = ncfile.findVariable("u");
            Variable vVar = ncfile.findVariable("v");
            Variable angVar = ncfile.findVariable("theta");

            // Get shape of theta variable:
            int[] angShape = angVar.getShape();
            
            // Get shape of u variable:
            int[] shape = uVar.getShape();
            for (int i=0; i<shape.length; i++)
                System.out.println("shape["+(i)+"]: "+shape[i]);

            ArrayDouble.D1 thetasD = (ArrayDouble.D1)angVar.read(new int[]{0}, angShape);


            thetas = new double[angShape[0]];
            for (int i=0; i<angShape[0]; i++)
                thetas[i] = thetasD.get(i);
            d_theta = thetas[1]-thetas[0];

            /*for (int i=0; i<angShape[0]; i++)
                System.out.println(thetas[i]);*/

            cfieldsU = new double[shape[1]][shape[2]][shape[0]];
            cfieldsV = new double[shape[1]][shape[2]][shape[0]];
            for (int k=0; k<shape[0]; k++) {
                ArrayDouble.D3 valueU = (ArrayDouble.D3)uVar.read(new int[] {k, 0, 0}, new int[] {1, shape[1], shape[2]});
                ArrayDouble.D3 valueV = (ArrayDouble.D3)vVar.read(new int[] {k, 0, 0}, new int[] {1, shape[1], shape[2]});
                for (int i = 0; i < shape[1]; i++) {
                    for (int j = 0; j < shape[2]; j++) {
                        cfieldsU[i][j][k] = valueU.get(0, i, j);
                        cfieldsV[i][j][k] = valueV.get(0, i, j);
                        //System.out.print(cfields[i][j][0] + " ");
                    }
                    //System.out.println("");
                }
            }



            ncfile.close();

            /*double angle;
            angle = 0; System.out.println("angle: "+angle+", index: "+chooseField(angle));
            angle = 5; System.out.println("angle: "+angle+", index: "+chooseField(angle));
            angle = 105; System.out.println("angle: "+angle+", index: "+chooseField(angle));
            angle = 10.5; System.out.println("angle: "+angle+", index: "+chooseField(angle));
            angle = -1; System.out.println("angle: "+angle+", index: "+chooseField(angle));
            angle = 359; System.out.println("angle: "+angle+", index: "+chooseField(angle));
            angle = -11; System.out.println("angle: "+angle+", index: "+chooseField(angle));

            double[][][][] test = new double[shape[1]][shape[2]][3][3];
            double[] directions = new double[] {135, 45, 270};
            double[] speeds = new double[] {5, 2, 3};
            setCurrentField(test, speeds, directions);

            File csvOutputFile = new File("cfield_u.csv");
            try (PrintWriter pw = new PrintWriter(csvOutputFile)) {
                for (int i=0; i<shape[1]; i++) {
                    for (int j = 0; j < shape[2]; j++) {
                        pw.print(test[i][j][0][0]+"\t");
                    }
                    pw.println("");
                }
                pw.close();
            }
            csvOutputFile = new File("cfield_v.csv");
            try (PrintWriter pw = new PrintWriter(csvOutputFile)) {
                for (int i=0; i<shape[1]; i++) {
                    for (int j = 0; j < shape[2]; j++) {
                        pw.print(test[i][j][0][1]+"\t");
                    }
                    pw.println("");
                }
                pw.close();
            }*/

        } catch (
                IOException ex) {
            ex.printStackTrace();

        } catch (
                InvalidRangeException e) {
            e.printStackTrace();

        }
    }

    /**
     * Set a full 3d current field based on directions and speeds given per layer. The directions are used to choose a
     * 2D field per layer, which is scaled according to the layer speed.
     * @param field
     * @param speeds
     * @param angles
     */
    public void setCurrentField(double[][][][] field, double[] speeds, double[] angles) {
        for (int k=0; k<angles.length; k++) {
            int fieldI = chooseField(angles[k]);
            for (int i=0; i<field.length; i++)
                for (int j=0; j<field[i].length; j++) {
                    //System.out.println("i="+i+", j="+j+", fieldI="+fieldI+", field.length="+field.length);
                    field[i][j][k][0] = speeds[k]*cfieldsU[i][j][fieldI];
                    field[i][j][k][1] = speeds[k]*cfieldsV[i][j][fieldI];
                    //field[i][j][k][0] = cfieldsU[i][j][fieldI];
                    //field[i][j][k][1] = cfieldsV[i][j][fieldI];
                    /*if (i==5 && j==15) {
                        System.out.println(field[i][j][k][0] + " / "+field[i][j][k][1]);
                        System.out.println();
                    }*/
                }
        }
    }

    /**
     * Given a theta angle, choose the index for the current field that best matches the angle.
     * @param theta The current angle (degrees)
     * @return The index of the best-matching field.
     */
    private int chooseField(double theta) {
        //System.out.println("theta (internal prior)="+theta);
        while (theta < 0)
            theta += 360;
        while (theta >= 360)
            theta -= 360;
        //System.out.println("theta (internal)="+theta);
        // Find closest value in thetas:
        double best = Double.MAX_VALUE;
        int bestI = -1;
        for (int i=0; i<thetas.length; i++) {
            if (Math.abs(theta-thetas[i]) < best) {
                bestI = i;
                best = Math.abs(theta-thetas[i]);
            }
        }
        return bestI;
    }

}
