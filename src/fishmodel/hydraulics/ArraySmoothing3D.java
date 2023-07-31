package fishmodel.hydraulics;

public class ArraySmoothing3D {

    // Function to apply Gaussian blur to a 3D array
    public static double[][][] applyGaussianBlur3D(double[][][] inputArray, int radius) {
        int width = inputArray.length;
        int height = inputArray[0].length;
        int depth = inputArray[0][0].length;
        double[][][] outputArray = new double[width][height][depth];

        // Generate a 3D Gaussian kernel
        double[][][] kernel = generateGaussianKernel3D(radius);

        // Apply Gaussian blur to each cell of the input array
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                for (int z = 0; z < depth; z++) {
                    double sum = 0.0;
                    double totalWeight = 0.0;

                    // Apply the kernel to neighboring cells
                    for (int i = -radius; i <= radius; i++) {
                        for (int j = -radius; j <= radius; j++) {
                            for (int k = -radius; k <= radius; k++) {
                                int newX = x + i;
                                int newY = y + j;
                                int newZ = z + k;

                                if (newX >= 0 && newX < width && newY >= 0 && newY < height && newZ >= 0 && newZ < depth) {
                                    double weight = kernel[i + radius][j + radius][k + radius];
                                    sum += inputArray[newX][newY][newZ] * weight;
                                    totalWeight += weight;
                                }
                            }
                        }
                    }

                    // Calculate the average value and set it to the output array
                    outputArray[x][y][z] = sum / totalWeight;
                }
            }
        }

        return outputArray;
    }

    // Function to generate a 3D Gaussian kernel with given radius
    private static double[][][] generateGaussianKernel3D(int radius) {
        int size = 2 * radius + 1;
        double sigma = radius / 3.0; // Adjust sigma for desired smoothing effect

        double[][][] kernel = new double[size][size][size];
        double twoSigmaSquare = 2 * sigma * sigma;
        double normalizingConstant = 0.0;

        for (int x = -radius; x <= radius; x++) {
            for (int y = -radius; y <= radius; y++) {
                for (int z = -radius; z <= radius; z++) {
                    kernel[x + radius][y + radius][z + radius] = Math.exp(-(x * x + y * y + z * z) / twoSigmaSquare);
                    normalizingConstant += kernel[x + radius][y + radius][z + radius];
                }
            }
        }

        // Normalize the kernel to ensure the sum of all elements is 1
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                for (int k = 0; k < size; k++) {
                    kernel[i][j][k] /= normalizingConstant;
                }
            }
        }

        return kernel;
    }

    // Sample usage
    public static void main(String[] args) {
        double[][][] inputArray = {
                {
                        {1.0, 2.0, 3.0},
                        {4.0, 5.0, 6.0},
                        {7.0, 8.0, 9.0}
                },
                {
                        {10.0, 11.0, 12.0},
                        {13.0, 14.0, 15.0},
                        {16.0, 17.0, 18.0}
                }
        };

        int radius = 1;
        double[][][] smoothedArray = applyGaussianBlur3D(inputArray, radius);

        // Print the smoothed array
        for (int i = 0; i < smoothedArray.length; i++) {
            for (int j = 0; j < smoothedArray[i].length; j++) {
                for (int k = 0; k < smoothedArray[i][j].length; k++) {
                    System.out.print(smoothedArray[i][j][k] + " ");
                }
                System.out.println();
            }
            System.out.println();
        }
    }
}
