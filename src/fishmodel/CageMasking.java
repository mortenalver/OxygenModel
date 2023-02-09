package fishmodel;

import java.util.ArrayList;

/**
 * Utility class for making a mask array. The mask array mirrors the tank grid, and contains true for all grid cells that are
 * inside the tank, and false for all cells outside. This can be used to prevent advection and diffusion through cage walls.
 */
public class CageMasking {

    /**
     * Create masking for a straight circular cylindrical tank.
     * @param dims Grid dimensions
      * @param dxy Horizontal grid size
      * @param radius Radius of cylinder
      * @return
     */
    public static boolean[][][] circularMasking(int[] dims, double dxy, double radius, boolean maskBottom) {
        boolean[][][] mask = new boolean[dims[0]][dims[1]][dims[2]];
        // Find center of grid:
        double[] center = new double[2];
        for (int i=0; i<2; i++)
            center[i] = ((double)(dims[i]-1))/2.0;

        for (int i=0; i<dims[0]; i++)
            for (int j=0; j<dims[1]; j++) {
                double distX = (double)i - center[0],
                            distY = (double)j- center[1];
                double rpos = dxy*Math.sqrt(distX*distX + distY*distY);
                boolean inside = rpos <= radius;
                for (int k=0; k<dims[2]-1; k++)
                    mask[i][j][k] = inside;
                if (maskBottom)
                    mask[i][j][dims[2]-1] = false;
                else mask[i][j][dims[2]-1] = inside;
            }

        return mask;
    }


    /**
     * Create masking for a cylindro-conical tank (typical fish cage shape).
     * @param dims Grid dimensions
     * @param dxy Horizontal grid size
     * @param dz Vertical grid size
     * @param radius Radius of cylindrical part
     * @param cylDepth Depth of cylindrical part
     * @param totDepth Depth of lowermost point of the cage
     * @return
     */
    public static boolean[][][] cylindroConicalMasking(int[] dims, double dxy, double dz, double radius, double cylDepth, double totDepth) {
        boolean[][][] mask = new boolean[dims[0]][dims[1]][dims[2]];
        // Find center of grid:
        double[] center = new double[2];
        for (int i=0; i<2; i++)
            center[i] = ((double)(dims[i]-1))/2.0;

        for (int i=0; i<dims[0]; i++)
            for (int j=0; j<dims[1]; j++) {
                double distX = (double)i - center[0],
                            distY = (double)j- center[1];
                double rpos = dxy*Math.sqrt(distX*distX + distY*distY);
                for (int k=0; k<dims[2]-1; k++) {
                    double depth = k*dz;
                    double rHere;
                    if (depth < cylDepth)
                        rHere = radius;
                    else if (depth < totDepth) {
                        rHere = radius*(totDepth-depth)/(totDepth-cylDepth);
                    }
                    else rHere = 0;
                    boolean inside = rHere > 0 && rpos <= rHere;
                    mask[i][j][k] = inside;

                }
            }

        return mask;
    }

    /**
     * Create masking for farm containing a grid of cages with given positions and radii
     * @param dims Grid dimensions
     * @param dxy Horizontal grid size
     * @param radius Radius of cylinder
     * @return
     */
    public static boolean[][][] fullFarmMasking(int[] dims, double dxy, ArrayList<double[]> cagePos, double radius, boolean maskBottom) {
        boolean[][][] mask = new boolean[dims[0]][dims[1]][dims[2]];
        // Initialize without cages first:
        for (int i=0; i<dims[0]; i++)
            for (int j=0; j<dims[1]; j++) {
                for (int k = 0; k < dims[2]; k++)
                    mask[i][j][k] = false;
            }

        for (double[] pos : cagePos) {
            for (int i = 0; i < dims[0]; i++)
                for (int j = 0; j < dims[1]; j++) {
                    double distX = (double) i - (pos[0]/dxy),
                            distY = (double) j - (pos[1]/dxy);
                    double rpos = dxy * Math.sqrt(distX * distX + distY * distY);
                    boolean inside = rpos <= radius;
                    if (inside) {
                        for (int k = 0; k < dims[2] - 1; k++)
                            mask[i][j][k] = true;
                    }
                }
        }

        return mask;
    }
}
