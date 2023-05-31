package fishmodel.enkf;

public class Util {
    public static double[] reshapeArray(double[][][] arr) {
        // Calculate the total number of elements in the 3D array
        int totalElements = arr.length * arr[0].length * arr[0][0].length;

        // Create a new 1D array with the same number of elements
        double[] reshapedArray = new double[totalElements];

        // Reshape the 3D array into a 1D array
        int index = 0;
        for (int i = 0; i < arr.length; i++) {
            for (int j = 0; j < arr[i].length; j++) {
                for (int k = 0; k < arr[i][j].length; k++) {
                    reshapedArray[index] = arr[i][j][k];
                    index++;
                }
            }
        }

        return reshapedArray;
    }

    public static int getStateIndex(int i, int j, int k, int[] cageDims) {
        return i*cageDims[1]*cageDims[2] + j*cageDims[2] + k;
    }
}
