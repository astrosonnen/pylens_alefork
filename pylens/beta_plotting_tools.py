import pylab
import numpy as np
import Image, ImageDraw, ImageFont
from yattaconfig import *


font = ImageFont.truetype("/usr/local/texlive/2015/texmf-dist/fonts/truetype/public/dejavu/DejaVuSans-Bold.ttf", 12)

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


def visual_comparison(data, models, sources=None):

    cuts = []
    resids = []
    i = 0
    for img in data:
        cut = np.percentile(img, 99.)
        cuts.append(cut)
        resids.append(img - models[i])
        i += 1

    #pylab.subplots_adjust(left=0., right=1., bottom=0., top=1., hspace=0., wspace=0.)
    pylab.subplot(2, 2, 1)
    pylab.imshow(make_rgbarray(data, cuts))
    pylab.xticks(())
    pylab.yticks(())

    pylab.subplot(2, 2, 2)
    pylab.imshow(make_rgbarray(models, cuts))
    pylab.xticks(())
    pylab.yticks(())

    pylab.subplot(2, 2, 3)
    pylab.imshow(make_rgbarray(resids, cuts))
    pylab.xticks(())
    pylab.yticks(())

    if sources is not None:
        pylab.subplot(2, 2, 4)
        pylab.imshow(make_rgbarray(sources, cuts))
        pylab.xticks(())
        pylab.yticks(())

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


def make_fail_curv_rgb(data, lenssub, arcimg, arccoords, arcmask, crapmask=None, fuzzmask=None, maskedge=None, \
                       outname='wrong_curvature.png'):

    cuts = []
    mask = []
    i = 0
    for img in data:
        cut = np.percentile(img, 99.)
        cuts.append(cut)
        if i==0:
            maskdata = arcmask*cut
            if crapmask is not None:
                maskdata += crapmask*cut
            if fuzzmask is not None:
                maskdata += 0.4*fuzzmask*cut
            mask.append(maskdata)
        else:
            mask.append(arcmask*cut)
        i += 1

    dlist = make_crazy_pil_format(data, cuts)
    lsublist = make_crazy_pil_format(lenssub, cuts)
    masklist = make_crazy_pil_format(mask, cuts)
    arclist = make_crazy_pil_format(arcimg, cuts)

    s = data[0].shape
    dim = Image.new('RGB', s, 'black')
    lsubim = Image.new('RGB', s, 'black')
    maskim = Image.new('RGB', s, 'black')
    arcim = Image.new('RGB', s, 'black')

    dim.putdata(dlist)
    lsubim.putdata(lsublist)
    maskim.putdata(masklist)
    arcim.putdata(arclist)

    x0 = s[1]/2
    y0 = s[0]/2
    if maskedge is not None:
        maskdraw = ImageDraw.Draw(maskim)
        maskdraw.ellipse((x0 - maskedge, y0 - maskedge, x0 + maskedge, y0 + maskedge), fill=None, outline='yellow')

    # draws the lines from the edges of the arc to the center

    def line_y(x, x1, x2, y1, y2):
        return (y2 - y1)/(x2 - x1)*(x - x2) + y2

    def line_x(y, x1, x2, y1, y2):
        return (x2 - x1)/(y2 - y1)*(y - y2) + x2

    def get_vertex(x1, x2, y1, y2):
        if x1 < 0.:
            xedge = 0.
            yedge = line_y(xedge, x1, x2, y1, y2)
            if yedge < 0.:
                yedge = 0.
                xedge = line_x(0., x1, x2, y1, y2)
            elif yedge > s[0]:
                yedge = s[0]
                xedge = line_x(yedge, x1, x2, y1, y2)
        elif x1 > s[1]:
            xedge = s[1]
            yedge = line_y(xedge, x1, x2, y1, y2)
            if yedge < 0.:
                yedge = 0.
                xedge = line_x(yedge, x1, x2, y1, y2)
            elif yedge > s[0]:
                yedge = s[0]
                xedge = line_x(yedge, x1, x2, y1, y2)
        else:
            if y1 < 0.:
                yedge = 0.
                xedge = line_x(yedge, x1, x2, y1, y2)
            elif y1 > s[0]:
                yedge = s[0]
                xedge = line_x(yedge, x1, x2, y1, y2)
            else:
                xedge = x1
                yedge = y1
        return (xedge, yedge)

    draw = ImageDraw.Draw(arcim)

    for coord in arccoords:
        xc, yc = coord[0]
        xa, ya = coord[1]
        xb, yb = coord[2]

        ca_coords = get_vertex(xc, xa, yc, ya)
        cb_coords = get_vertex(xc, xb, yc, yb)


        draw.line((ca_coords[0], s[1] - ca_coords[1], xa, s[1] - ya), fill=(255, 255, 255))
        draw.line((cb_coords[0], s[1] - cb_coords[1], xb, s[1] - yb), fill=(255, 255, 255))

    im = Image.new('RGB', (4*data[0].shape[0], data[0].shape[1]), 'black')

    im.paste(dim, (0, 0,))
    im.paste(lsubim, (s[1], 0))
    im.paste(maskim, (2*s[1], 0))
    im.paste(arcim, (3*s[1], 0))

    im.save(outname)


def make_full_rgb(candidate, image_set, maskedge=None, outname='full_model.png', nsig_cut=5., success=None):

    cuts = []
    rescuts = []
    data = []
    lenssub = []
    ringresid = []
    lensresid = []
    sersicresid = []
    lensmodel = []
    ringmodel = []
    sersicmodel = []
    source = []

    mask = []
    i = 0

    ncol = 6

    for band in rgbbands:
        img = candidate.sci[band]
        data.append(img)
        maskimg = 0.*img
        cut = np.percentile(img[candidate.R < 30.], 99.)
        cuts.append(cut)
        rescuts.append(np.percentile(candidate.lenssub_resid[band][candidate.R < 30.], 99.))
        lenssub.append(candidate.lenssub_resid[band])

        lmodel = 0.*img
        for mimg in candidate.lensfit_model[band]:
            lmodel += mimg
        lensmodel.append(lmodel)

        source.append(candidate.lensfit_model[band][-1])

        lensresid.append(candidate.sci[band] - lmodel)

        for image in image_set['images'] + image_set['foregrounds'] + image_set['bad_arcs']:
            maskimg[image['footprint'] > 0] = cut

        if i==0:
            for junk in image_set['junk']:
                maskimg[junk['footprint'] > 0] = cut
        elif i==1:
            for arc in image_set['arcs']:
                maskimg[arc['footprint'] > 0] = cut

        mask.append(maskimg)
        i += 1

    dlist = make_crazy_pil_format(data, cuts)
    lsublist = make_crazy_pil_format(lenssub, rescuts)
    masklist = make_crazy_pil_format(mask, cuts)
    slist = make_crazy_pil_format(source, cuts)

    lmlist = make_crazy_pil_format(lensmodel, cuts)
    lrlist = make_crazy_pil_format(lensresid, cuts)

    s = (data[0].shape[1], data[0].shape[0])
    dim = Image.new('RGB', s, 'black')
    lsubim = Image.new('RGB', s, 'black')
    maskim = Image.new('RGB', s, 'black')
    sim = Image.new('RGB', s, 'black')
    lmim = Image.new('RGB', s, 'black')
    lrim = Image.new('RGB', s, 'black')

    dim.putdata(dlist)
    lsubim.putdata(lsublist)
    maskim.putdata(masklist)
    sim.putdata(slist)
    lmim.putdata(lmlist)
    lrim.putdata(lrlist)

    if len(candidate.ringfit_model) > 0:
        rcol = ncol
        ncol += 2
        for band in rgbbands:
            rmodel = 0.*img
            for mimg in candidate.ringfit_model[band]:
                rmodel += mimg
            ringmodel.append(rmodel)
            ringresid.append(candidate.sci[band] - rmodel)

        rmlist = make_crazy_pil_format(ringmodel, cuts)
        rrlist = make_crazy_pil_format(ringresid, cuts)

        rmim = Image.new('RGB', s, 'black')
        rrim = Image.new('RGB', s, 'black')

        rmim.putdata(rmlist)
        rrim.putdata(rrlist)

    if len(candidate.sersicfit_model) > 0:
        scol = ncol
        ncol += 2
        for band in rgbbands:
            smodel = 0.*img
            for mimg in candidate.sersicfit_model[band]:
                smodel += mimg
            sersicmodel.append(smodel)
            sersicresid.append(candidate.sci[band] - smodel)

        cmlist = make_crazy_pil_format(sersicmodel, cuts)
        crlist = make_crazy_pil_format(sersicresid, cuts)

        cmim = Image.new('RGB', s, 'black')
        crim = Image.new('RGB', s, 'black')

        cmim.putdata(cmlist)
        crim.putdata(crlist)

    x0 = s[1]/2
    y0 = s[0]/2
    if maskedge is not None:
        maskdraw = ImageDraw.Draw(maskim)
        maskdraw.ellipse((x0 - maskedge, y0 - maskedge, x0 + maskedge, y0 + maskedge), fill=None, outline='yellow')

    im = Image.new('RGB', (ncol*data[0].shape[0], data[0].shape[1]), 'black')

    im.paste(dim, (0, 0,))
    im.paste(lsubim, (s[1], 0))
    im.paste(maskim, (2*s[1], 0))
    im.paste(lmim, (3*s[1], 0))
    im.paste(sim, (4*s[1], 0))
    im.paste(lrim, (5*s[1], 0))

    draw = ImageDraw.Draw(im)
    draw.text((10, s[0] - 20), 'HSCJ'+candidate.name, font=font, fill='white')
    draw.text((10 + 5*s[1], s[0] - 20), '%2.1f'%candidate.lensfit_chi2, font=font, fill='white')

    if len(candidate.ringfit_model) > 0:
        im.paste(rmim, (rcol*s[1], 0))
        im.paste(rrim, ((rcol+1)*s[1], 0))
        draw.text((10 + rcol*s[1], s[0] - 20), '%2.1f'%candidate.ringfit_chi2, font=font, fill='white')

    if len(candidate.sersicfit_model) > 0:
        im.paste(cmim, (scol*s[1], 0))
        im.paste(crim, ((scol+1)*s[1], 0))
        draw.text((10 + scol*s[1], s[0] - 20), '%2.1f'%candidate.sersicfit_chi2, font=font, fill='white')

    if success is not None:
        if success:
            draw.ellipse((10, 10, 30, 30), fill=None, outline=(0, 255, 0))
            draw.ellipse((11, 11, 29, 29), fill=None, outline=(0, 255, 0))
        else:
            draw.line((10, 10, 30, 30), fill=(255, 0, 0), width=2)
            draw.line((30, 10, 10, 30), fill=(255, 0, 0), width=2)

    im.save(outname)
