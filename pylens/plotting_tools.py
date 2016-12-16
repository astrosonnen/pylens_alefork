import pylab
import numpy as np
import Image, ImageDraw, ImageFont
from yattaconfig import *



def make_rgbarray(images, cuts):

    scaled = []
    for i in range(3):
        img = images[i].copy()

        img[img<0.] = 0.
        img *= 255./cuts[i]
        img[img>255.] = 255.
        img = np.uint8(img.round())
        img = np.flipud(img)
        scaled.append(img.T)

    rgbarray = np.array(scaled).T
    return rgbarray


def make_crazy_pil_format(data, cuts):

    newlist = []
    for i in range(0, 3):
        flatd = np.flipud(data[i]).flatten()
        flatd[flatd<0.] = 0.
        flatd *= 255./cuts[i]
        flatd[flatd>255.] = 255.
        flatd = np.uint8(flatd.round())
        newlist.append(flatd)

    l = []

    for i in range(0, data[0].size):
        l.append((newlist[0][i], newlist[1][i], newlist[2][i]))

    return l


def make_model_rgb(sci, light_model, source_model, cuts=(99., 99., 99.), outname='model_rgb.png'):

    auto_cuts = []
    data = []
    lensresid = []
    lensmodel = []

    i = 0

    ncol = 4

    for i in range(3):
        data.append(sci[i])
        cut = np.percentile(sci[i], cuts[i])
        auto_cuts.append(cut)

        lensmodel.append(light_model[i] + source_model[i])

        lensresid.append(sci[i] - light_model[i] - source_model[i])

    dlist = make_crazy_pil_format(data, auto_cuts)
    slist = make_crazy_pil_format(source_model, auto_cuts)

    lmlist = make_crazy_pil_format(lensmodel, auto_cuts)
    lrlist = make_crazy_pil_format(lensresid, auto_cuts)

    s = (data[0].shape[1], data[0].shape[0])
    dim = Image.new('RGB', s, 'black')
    sim = Image.new('RGB', s, 'black')
    lmim = Image.new('RGB', s, 'black')
    lrim = Image.new('RGB', s, 'black')

    dim.putdata(dlist)
    sim.putdata(slist)
    lmim.putdata(lmlist)
    lrim.putdata(lrlist)

    im = Image.new('RGB', (ncol*data[0].shape[0], data[0].shape[1]), 'black')

    im.paste(dim, (0, 0,))
    im.paste(lmim, (1*s[1], 0))
    im.paste(sim, (2*s[1], 0))
    im.paste(lrim, (3*s[1], 0))

    im.save(outname)

