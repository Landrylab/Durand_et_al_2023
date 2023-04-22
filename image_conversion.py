# Conversion of images taken by S & P robot

def get_modified_image(image, outp, pic_fname, xoff=500, yoff=300):
    
    from pathlib import Path
    import os, os.path
    from skimage import io,util
   
    # xoff/yoff = number of pixels that will be cropped from the picture
    
    # Output folder is created only if it doesn't already exist
    Path(outp).mkdir(parents=True, exist_ok=True)
    
    # Picture is cropped
    ## I've added +80 because the plate is not perfectly centered, which means there's more to crop on the left than on the right
    y1,y2,x1,x2=yoff,image.shape[0]-yoff,xoff,image.shape[1]-xoff+80
    cropped = image[y1:y2,x1:x2]

    # Picture is converted into levels of gray
    from skimage.color import rgb2gray
    grayscale = rgb2gray(cropped)
    
    # Levels of gray are inverted: gray on white background equals growth
    grayscale=util.invert(grayscale)
    out=util.img_as_ubyte(grayscale)

    # Converted picture is saved
    io.imsave(outp+pic_fname, out)
    
    return out, outp+pic_fname

# Batch conversion
for f in os.listdir(p):
    fname = os.path.basename(f)
    image = io.imread(p+fname)
    get_modified_image(image, outp, pic_fname = fname)