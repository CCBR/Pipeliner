library("magick")

tree <- image_read("sample_network.bmp")
image_write(tree, path = "sample_network_mqc.png", format = "png")