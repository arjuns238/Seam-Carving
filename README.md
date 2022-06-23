# Seam-Carving

Seam-carving is a content-aware image resizing technique where the image is reduced in size by one pixel of height (or width) at a time. A vertical seam in an image is a path of pixels connected from the top to the bottom with one pixel in each row; a horizontal seam is a path of pixels connected from the left to the right with one pixel in each column.

I implemented this resizing technique using Dijkstra's shortest path algorithm. The pixels in the image were mapped to an energy function which thereby gave the weights of the graph. Each pixel was connected to the 3 closest pixels in the next row. Using the shortest path algorithm, I found the lowest energy pixels that connected the top row to the bottom row and removed these pixels from the image.
