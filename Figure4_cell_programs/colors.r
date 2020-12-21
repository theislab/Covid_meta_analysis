require(RColorBrewer)

desat = function(cols, sat=0.5) {
    X = diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
    hsv(X[1,], X[2,], X[3,])
}

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

set.colors = c(brewer.pal(9, 'Set1'), brewer.pal(7, 'Set2'), brewer.pal(12, 'Set3')[c(3,9,8)], 'violetred4')
set.colors[6] = 'khaki2'
set.colors[8] = 'lightskyblue2'
set.colors = rep(set.colors, 10)

tsne.colors = c(brewer.pal(7, 'Set2'), brewer.pal(9, 'Set1'), brewer.pal(12, 'Set3')[c(3,9,8)], 'violetred4')
tsne.colors[13] = 'khaki2'
tsne.colors[15] = 'lightskyblue2'
tsne.colors = rep(tsne.colors, 10)
tsne.colors = desat(tsne.colors, sat=.75)

material.heat <- function(n)
{
    mh = c(
        #"#607D8B", #blue grey
        "#283593", #indigo 800
        "#3F51B5",  #indigo
        "#2196F3", # blue
        #"#03A9F4", # light blue
        "#00BCD4", #cyan
        #"#009688", # teal
        "#4CAF50", #green
        "#8BC34A", #light green
        "#CDDC39", # lime
        "#FFEB3B", #yellow
        "#FFC107", # amber
        "#FF9800", # organe
        "#FF5722", # deep orange)
        "#f44336")
    colorRampPalette(mh)(n)
}

nmf.colors = c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")

tol1qualitative=c("#4477AA")
tol2qualitative=c("#4477AA", "#CC6677")
tol3qualitative=c("#4477AA", "#DDCC77", "#CC6677")
tol4qualitative=c("#4477AA", "#117733", "#DDCC77", "#CC6677")
tol5qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677")
tol6qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499")
tol7qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499")
tol8qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677","#AA4499")
tol9qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")
tol10qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
tol11qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
tol12qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499")


# material colors
# ---------------

material = read.table('~/code/single_cell/material_colors.txt', row.names=1, comment.char='', stringsAsFactors=F)
rownames(material) = tolower(rownames(material))

shades = function(color, n){
    colorRampPalette(material[color,])(n)
}
disc.colors = function(n, shade=6){
    disc.order = c('red', 'light_blue', 'amber', 'green', 'grey')
    if(n <= length(disc.order)){
        material[disc.order[1:n], shade]
    } else if(n <= nrow(material)){
        material[as.integer(seq(from=1, to=nrow(material), length.out=n)), shade]
    } else {
        colorRampPalette(material[,shade])(n)
    }
}
cont.colors = function(shade=6){
    colorRampPalette(material[,shade])(100)
}
name.colors = function(names, shade=6){
    material[names, shade]
}

set3.colors = c('#80b1d3', '#fdb462', '#b3de69')
set3.colors = c(set3.colors, setdiff(brewer.pal(12, 'Set3'), set3.colors))


# get similar colors
# ------------------

col2hcl <- function(x, maxColorValue=255, ...){
    library(colorspace)
    x = RGB(t(col2rgb(x))[,1:3,drop=F]/maxColorValue)
    x = colorspace::coords(as(x, 'polarLUV'))
    x[is.na(x)] = 0
    x[1,]
}

nice_colors = function(n, col=NULL, type='fancy', hmin=NULL, hmax=NULL, cmin=NULL, cmax=NULL, lmin=NULL, lmax=NULL, plot=F){
    
    library(grDevices)
    library(hues)
    
    # initialize palette
    beg = end = c()

    # map color to hues
    if(!is.null(col)){
        h0 = col2hcl(col)[['H']]
	hmin = h0 - 30
	hmax = h0 + 30
    }

    # get preset color schemes
    if(!is.null(type)){
	if(type == 'single'){cmin = 10; cmax = 100; lmin = 35; lmax = 100}	
        if(type == 'shades'){hmin = 0; hmax = 240; cmin = 0; cmax = 15; lmin = 0; lmax = 100;}
	if(type == 'tarnish'){hmin = 0; hmax = 360; cmin = 0; cmax = 15; lmin = 30; lmax = 70;}
	if(type == 'pastel'){hmin = 0; hmax = 360; cmin = 0; cmax = 30; lmin = 70; lmax = 100;}
	if(type == 'pimp1'){hmin = 0; hmax = 360; cmin = 30; cmax = 100; lmin = 25; lmax = 70;}
	if(type == 'pimp2'){hmin = 0; hmax = 360; cmin = 30; cmax = 100; lmin = 55; lmax = 100;}
	if(type == 'intense'){hmin = 0; hmax = 360; cmin = 20; cmax = 100; lmin = 15; lmax = 80;}
	if(type == 'fluo'){hmin = 0; hmax = 300; cmin = 35; cmax = 100; lmin = 75; lmax = 100;}
	if(type == 'fancy'){hmin = 0; hmax = 360; cmin = 15; cmax = 40; lmin = 70; lmax = 100;}
	if(type == 'large'){hmin = 0; hmax = 360; cmin = 0; cmax = 85; lmin = 30; lmax = 95;}
    }
    
    # fix chroma values
    cmin = cmin*180/100.
    cmax = cmax*180/100.
    
    # handle edge cases
    if((hmin < 0 & hmax < 0) | (hmin > 360 & hmax > 360)){stop('error: invalid hmin and hmax range')}
    if(hmin < 0){
        n1 = round(n*(-hmin)/(hmax - hmin))
	if(n1 > 0){
	    beg = iwanthue(n1, hmin=360+hmin, hmax=360, cmin=cmin, cmax=cmax, lmin=lmin, lmax=lmax)
	    n = n - n1	    
	}
	hmin = 0
    }
    if(hmax > 360){
        n1 = round(n*(hmax - 360)/(hmax - hmin))
	if(n1 > 0){
	    end = iwanthue(n1, hmin=0, hmax=hmax-360, cmin=cmin, cmax=cmax, lmin=lmin, lmax=lmax)
	    n = n - n1
	}
	hmax = 360
    }
    pal = c(beg, iwanthue(n, hmin=hmin, hmax=hmax, cmin=cmin, cmax=cmax, lmin=lmin, lmax=lmax), end)
    unname(pal)
}

set2.colors = c(brewer.pal(9, 'Set1'), brewer.pal(7, 'Set2'), brewer.pal(12, 'Set3')[c(3,9,8)], 'violetred4')
set2.colors[6] = 'khaki2'
set2.colors[8] = 'lightskyblue2'
set2.colors = c(set2.colors, nice_colors(16, type='intense'), nice_colors(16, type='large'))

