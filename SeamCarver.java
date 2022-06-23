/* *****************************************************************************
 *  Name:
 *  Date:
 *  Description:
 **************************************************************************** */

import edu.princeton.cs.algs4.Picture;

public class SeamCarver {
    // create a seam carver object based on the given picture
    private Picture pic;
    private double[][] energy;
    private double[][] distTo;
    private boolean isTransposed;
    private int[][] pixels;
    private boolean isTransposedPixels;


    public SeamCarver(Picture picture) {
        if (picture == null)
            throw new IllegalArgumentException("Picture is null");
        pic = new Picture(picture);
        energy = new double[picture.height()][picture.width()];
        isTransposed = false;
        isTransposedPixels = false;
        pixels = new int[picture.height()][picture.width()];
        // In image processing, pixel(i, j) refers to the pixel in column i and row j
        for (int row = 0; row < picture.height(); row++) {
            for (int col = 0; col < picture.width(); col++) {
                pixels[row][col] = pic.getRGB(col, row);
            }
        }
        for (int row = 0; row < picture.height(); row++) {
            for (int col = 0; col < picture.width(); col++) {
                energy[row][col] = CalculateEnergy(col, row);
            }
        }
    }


    private double CalculateEnergy(int col, int row) {
        // if pixel is an edge case then weight is 1000
        if (col == 0 || row == 0 || col == pixels[0].length - 1 || row == pixels.length - 1)
            return 1000;

        double rDiff = Math.abs(getRedval(col + 1, row) - getRedval(col - 1, row));
        double bDiff = Math.abs(getBlueval(col + 1, row) - getBlueval(col - 1, row));
        double gDiff = Math.abs(getGreenval(col + 1, row) - getGreenval(col - 1, row));
        double ryDiff = Math.abs(getRedval(col, row + 1) - getRedval(col, row - 1));
        double byDiff = Math.abs(getBlueval(col, row + 1) - getBlueval(col, row - 1));
        double gyDiff = Math.abs(getGreenval(col, row + 1) - getGreenval(col, row - 1));

        double xGrad = Math.pow((rDiff), 2) + Math.pow(bDiff, 2) + Math.pow(gDiff, 2);
        double yGrad = Math.pow((ryDiff), 2) + Math.pow(byDiff, 2) + Math.pow(gyDiff, 2);


        return Math.sqrt(xGrad + yGrad);
    }


    private double getRedval(int col, int row) {
        int rgb = pixels[row][col];
        return (rgb >> 16) & 0xFF;
    }

    private double getBlueval(int col, int row) {
        int rgb = pixels[row][col];
        return (rgb >> 8) & 0xFF;
    }

    private double getGreenval(int col, int row) {
        int rgb = pixels[row][col];
        return (rgb >> 0) & 0xFF;
    }


    private boolean validateVertices(int col, int row) {
        return (col >= 0 && col < energy[0].length && row >= 0 && row < energy.length);

    }


    // current picture
    public Picture picture() {
        if (isTransposedPixels)
            transposePixels();
        isTransposedPixels = false;
        pic = new Picture(pixels[0].length, pixels.length); // check thissss
        for (int row = 0; row < pixels.length; row++) {
            for (int col = 0; col < pixels[0].length; col++) {
                pic.setRGB(col, row, pixels[row][col]);
            }
        }
        return pic;
    }

    // width of current picture
    public int width() {
        return energy[0].length;

    }

    // height of current picture
    public int height() {
        return energy.length;
    }

    // energy of pixel at column x and row y
    public double energy(int x, int y) {
        if (x < 0 || x >= energy[0].length || y < 0 || y >= energy.length)
            throw new IllegalArgumentException("invalid energy");
        return energy[y][x];
    }


    // sequence of indices for horizontal seam
    public int[] findHorizontalSeam() {
        if (!isTransposed)
            transpose();
        isTransposed = true;
        // energy.length = picture.height
        // energy[0].length = picture.width

        distTo = new double[energy.length][energy[0].length];


        for (int row = 0; row < energy.length; row++) {
            for (int col = 0; col < energy[0].length; col++) {
                if (row == 0)
                    distTo[row][col] = 0;       // experimental,
                else
                    distTo[row][col] = Double.POSITIVE_INFINITY;

            }
        }
        int[] horizontalSeam = findSeam();
        if (isTransposed)
            transpose();
        isTransposed = false;
        return horizontalSeam;
    }

    // sequence of indices for vertical seam
    public int[] findVerticalSeam() {


        if (isTransposed)
            transpose();
        isTransposed = false;
        distTo = new double[energy.length][energy[0].length];


        for (int row = 0; row < energy.length; row++) {
            for (int col = 0; col < energy[0].length; col++) {
                if (row == 0)
                    distTo[row][col] = 0;
                else
                    distTo[row][col] = Double.POSITIVE_INFINITY;

            }
        }


        return findSeam();

    }

    private int[] findSeam() {


        for (int row = 0; row < energy.length; row++) {   // height
            for (int col = 0; col < energy[0].length; col++) {  // width
                relax(col, row, col - 1, row + 1);
                relax(col, row, col, row + 1);
                relax(col, row, col + 1, row + 1);
            }
        }

        int minAtLastRow = minAtLastRow(distTo);
        return returnOrder(distTo, minAtLastRow);


    }


    private void relax(int col, int row, int childCol, int childRow) {
        // energy.length = picture.height
        // energy[0].length = picture.width
        if (validateVertices(childCol, childRow)) {
            if (distTo[childRow][childCol] > distTo[row][col] + energy[childRow][childCol])
                distTo[childRow][childCol] = distTo[row][col] + energy[childRow][childCol];
        }
    }

    private int minAtLastRow(double[][] arr) {
        int row = arr.length - 1;
        double min = Double.POSITIVE_INFINITY;
        int minIndex = -1;
        for (int col = 0; col < arr[0].length; col++) {
            if (arr[row][col] < min) {
                min = arr[row][col];
                minIndex = col;
            }
        }
        return minIndex;
    }

    private int[] returnOrder(double[][] arr, int minIndex) {

        boolean parent1;
        boolean parent2;
        boolean parent3;
        int[] order = new int[arr.length];
        for (int row = arr.length - 1; row >= 0; row--) {

            order[row] = minIndex;

            parent1 = validateVertices(minIndex - 1, row - 1);
            parent2 = validateVertices(minIndex, row - 1);
            parent3 = validateVertices(minIndex + 1, row - 1);
            if (parent1 && parent2 && parent3)
                minIndex = minimumIndex(minIndex - 1, minIndex, minIndex + 1, row - 1);
            else if (parent1 && parent2)
                minIndex = minimumIndex(minIndex - 1, minIndex, row - 1);
            else if (parent2 && parent3)
                minIndex = minimumIndex(minIndex, minIndex + 1, row - 1);
            else
                break;

        }


        return order;

    }

    private int minimumIndex(int col1, int col2, int col3, int row) {
        double val1 = distTo[row][col1];
        double val2 = distTo[row][col2];
        double val3 = distTo[row][col3];
        double min = Math.min((Math.min(val1, val2)), val3);
        if (min == val1)
            return col1;
        if (min == val2)
            return col2;
        else
            return col3;
    }

    private int minimumIndex(int col1, int col2, int row) {
        double val1 = distTo[row][col1];
        double val2 = distTo[row][col2];
        double min = Math.min(val1, val2);
        if (min == val1)
            return col1;
        else
            return col2;
    }


    private void transpose() {

        double[][] temp = new double[energy[0].length][energy.length];
        for (int i = 0; i < energy.length; i++) { // height
            for (int j = 0; j < energy[0].length; j++) {
                temp[j][i] = energy[i][j];
            }
        }
        energy = new double[temp.length][temp[0].length];
        for (int i = 0; i < temp[0].length; i++) { // height
            for (int j = 0; j < temp.length; j++) { // width
                energy[j][i] = temp[j][i];
            }

        }
    }

    private void transposePixels() {
        int[][] temp = new int[pixels[0].length][pixels.length];
        for (int i = 0; i < pixels.length; i++) { // height
            for (int j = 0; j < pixels[0].length; j++) {
                temp[j][i] = pixels[i][j];
            }
        }
        int row = pixels.length;
        int col = pixels[0].length;
        pixels = new int[col][row];
        for (int i = 0; i < temp[0].length; i++) { // height
            for (int j = 0; j < temp.length; j++) { // width
                pixels[j][i] = temp[j][i];
            }

        }
        temp = null;
    }


    // remove horizontal seam from current picture
    public void removeHorizontalSeam(int[] seam) {
        if (seam == null)
            throw new IllegalArgumentException("seam is null");
        if (!isTransposed) {
            transpose();
            isTransposed = true;
        }

        if (!isTransposedPixels) {
            transposePixels();
            isTransposedPixels = true;
        }
        removeVerticalSeam(seam);
        transpose();
        transposePixels();
        isTransposed = false;
        isTransposedPixels = false;


    }

    // remove vertical seam from current picture
    public void removeVerticalSeam(int[] seam) {
        if (seam == null)
            throw new IllegalArgumentException("seam is null");
        int current = seam[0];
        if (seam.length != energy.length)
            throw new IllegalArgumentException("seam argument invalid");
        for (int i = 0; i < seam.length; i++) {
            if (seam[i] < 0 || seam[i] >= energy[0].length)
                throw new IllegalArgumentException("seam argument invalid");
            if (Math.abs(seam[i] - current) > 1)
                throw new IllegalArgumentException("seam argument invalide");
            current = seam[i];
        }
        int[][] tempPixels = new int[energy.length][energy[0].length - 1];

        for (int row = 0; row < energy.length; row++) {
            for (int col = 0; col < energy[0].length; col++) {
                if (seam[row] < 0)
                    throw new IllegalArgumentException("seam argument invalid");
                // changing pic
                if (col == tempPixels[0].length - 1 && col == seam[row])
                    continue;
                if (col < seam[row])
                    tempPixels[row][col] = pixels[row][col];
                if (col > seam[row])
                    tempPixels[row][col - 1] = pixels[row][col];

            }
        }
        pixels = new int[tempPixels.length][tempPixels[0].length];
        for (int row = 0; row < tempPixels.length; row++) {
            for (int col = 0; col < tempPixels[0].length; col++) {
                pixels[row][col] = tempPixels[row][col];
            }
        }


        double[][] newEnergy = new double[energy.length][energy[0].length - 1];

        for (int row = 0; row < energy.length; row++) {
            for (int col = 0; col < energy[0].length; col++) {

                if (col < seam[row] - 1 && col < newEnergy[0].length)
                    newEnergy[row][col] = energy[row][col];
                if (col > seam[row] && col < newEnergy[0].length)
                    newEnergy[row][col] = energy[row][col + 1];
                if ((col == seam[row] - 1 || col == seam[row]) && col < newEnergy[0].length) {
                    newEnergy[row][col] = CalculateEnergy(col, row);
                }

            }
        }
        energy = new double[newEnergy.length][newEnergy[0].length];
        for (int row = 0; row < newEnergy.length; row++) {
            for (int col = 0; col < newEnergy[0].length; col++) {
                energy[row][col] = newEnergy[row][col];
            }
        }

    }

    public static void main(String[] args) {
    }
}
