#seurat()

make_edges = function(data, pairs, subset=NULL, binary=FALSE, symmetric=TRUE, max_cells_per_gene=2, verbose=TRUE){
    
    # this function creates a (cells x cells) edge matrix
    # weights wij are calculated as the product of the node weights, wij=wi*wj
    # input arguments:
    # data = (cells x genes) weights matrix or (cells, genes, weights) list
    # pairs = (gene, gene) list
    # subset = (gene -> cells) list

    library(Matrix)
    library(tidyr)
    library(data.table)
    
    # format data
    if(ncol(data) == 3){
        colnames(data) = c('cell', 'gene', 'value')
	data = data.frame(spread(data, gene, value, fill=0), row.names=1)
	data = as(as.matrix(data), 'sparseMatrix')
    }
    print(dim(data))
    
    # make binary
    if(binary == FALSE & max(data) > 1){
        print('[warning] max(data) > 1')
    }
    
    if(binary == TRUE){
        data = data > 0
    }

    # remove promiscuous genes
    data = data[,colSums(data > 0, na.rm=T) <= max_cells_per_gene]
    print(dim(data))
    
    # make adjacency matrix
    cells = rownames(data)
    A = matrix(0, nrow=length(cells), ncol=length(cells))
    rownames(A) = colnames(A) = cells
    
    # make labels
    L = sapply(cells, function(a) sapply(cells, function(b) c(), simplify=F), simplify=F)
    
    # add edges
    count = 0
    for(i in 1:nrow(pairs)){
        u = as.character(pairs[i,1]) # source gene
	v = as.character(pairs[i,2]) # target gene
	if(u %in% colnames(data) & v %in% colnames(data)){
	    
	    Ai = data[,u] %*% t(data[,v])
	    
	    if(!is.null(subset)){

	        # subset on gene u
		Au = Ai
		i = rownames(data) %in% subset[[u]]
		Au[!i,] = 0
		
		# subset on gene v
		Av = Ai
		i = rownames(data) %in% subset[[v]]
		Av[!i,] = 0
		
		# sum components
		Ai = Au + Av

		# skip on empty
		if(sum(Ai) == 0){
		    next
		}
	    }
	    
	    if(symmetric == TRUE){
	        Ai = Ai + t(Ai)
	    }
	    
	    A = A + Ai
	    
	    # update labels
	    cells1 = rownames(data)[data[,u] > 0]
	    cells2 = rownames(data)[data[,v] > 0]
	    for(ca in cells1){
	        for(cb in cells2){
		    L[[ca]][[cb]] = c(L[[ca]][[cb]], paste(u, v, sep='-'))
		    if(symmetric == TRUE){
		        L[[cb]][[ca]] = c(L[[cb]][[ca]], paste(u, v, sep='-'))
		    }
		}
	    }
	    
	    # print message
	    if(verbose == TRUE){
	        print(paste0('Edge: ', u, ', ', v, ' (', paste(cells1, collapse=', '), ') <--> (', paste(cells2, collapse=', '), ')'))
	    }
	    count = count + 1
        }
    }
    
    # fix labels
    for(ca in names(L)){
        for(cb in names(L[[ca]])){
	    L[[ca]][[cb]] = unique(L[[ca]][[cb]])
	}
    }
    
    print(paste(count, 'total interactions'))
    return(list(edges=A, labels=L))
}


rescale_vector = function(x, target=c(0, 1), f=identity, abs=FALSE){
    # map vector onto target interval, then apply f
    if(min(target) == max(target)){
        return(rep(max(target), length(x)))
    }
    a = target[[1]]
    b = target[[2]]
    if(min(x) == max(x)){
        rep(max(target), length(x))
    } else {
        if(abs == TRUE){
            x = (x - min(abs(x)))/(max(abs(x)) - min(abs(x)))*(b - a) + a
	} else {
            x = (x - min(x))/(max(x) - min(x))*(b - a) + a
	}
        f(x)
    }
}

plot_network = function(mat=NULL, edges=NULL, lab=NULL, f_weights=identity, rescale_weights=c(0,1), f_coords=identity, rescale_coords=c(0,1), alpha=1, remove_self=FALSE,
                        edge.label.cex=1, vertex.label.cex=1, node_color=NULL){
    
    # plot network diagram from matrix or edgelist
    # mat = (source x target) weights matrix
    # edges = (source, target, weight) edgelist
    # lab = list(source -> target -> label)
    
    library(igraph)
    library(data.table)
    library(tidyverse)
    
    # convert matrix to edges
    if(!is.null(mat)){
        edges = gather(as.data.frame(mat) %>% rownames_to_column('source'), target, weight, -source)
    }
    
    # select edges
    edges = data.frame(as.data.table(edges)[,.SD[which.max(abs(weight))],.(source, target)], stringsAsFactors=F)
    edges$source = as.character(edges$source)
    edges$target = as.character(edges$target)
    
    # remove zeros and de-duplicate
    edges = edges[edges$weight != 0,]
    i = apply(edges, 1, function(a) paste(sort(c(a[[1]], a[[2]])), collapse=' '))
    edges = edges[!duplicated(i),]
    
    # modify edges
    edges$weight = f_weights(edges$weight)
    
    # delete zero and self edges
    edges = edges[abs(edges$weight) > 0,]
    if(remove_self == TRUE){
        edges = edges[edges$source != edges$target,]
    }
    
    # get nodes
    nodes = data.frame(id=sort(unique(c(edges[,1], edges[,2]))), stringsAsFactors=F)
    
    # adjust edges
    if(all(edges$weight > 0)){
        edges$color = '#000000'
    } else {
        edges$color = ifelse(edges$weight > 0, 'blue', 'red')
        edges$weight = abs(edges$weight)	
    }
    
    # make graph
    g = graph.data.frame(d=edges, vertices=nodes, directed=FALSE)
    if(!is.null(node_color)){
        V(g)$color = node_color[names(V(g))]
    } else {
        V(g)$color = '#cccccc'
    }
    l = layout_with_fr(g, weight=abs(edges$weight) > 0)
    
    # make labels
    lab = cbind(l, nodes)
    colnames(lab) = c('x', 'y', 'lab')
    
    # rescale coordinates
    l = apply(l, 2, function(a) rescale_vector(a, target=rescale_coords, f=f_coords, abs=TRUE))
    
    # rescale edge widths
    w = edges$weight
    w = rescale_vector(w, target=rescale_weights, f=f_weights, abs=TRUE)    
    
    # plot and save
    plot(g, edge.width=w, vertex.color=V(g)$color, edge.color=E(g)$color, edge.label=E(g)$label, edge.label.color='black', vertex.label.color='black', layout=l,
         edge.label.cex=edge.label.cex, vertex.label.cex=vertex.label.cex)

}


ggplot_network = function(graph=NULL, edges=NULL, node_sizes=NULL, node_colors=NULL, edge_color='#cccccc', edge_colors=NULL, scale_nodes=c(4,4), scale_edges=c(.5,.5), curvature=0,
                          smin=NULL, smax=NULL,
                          alpha=1, symm=TRUE, do.legend=TRUE, legend_title=NULL, out=NULL, nrow=1, ncol=1, qmax=1, ggtitle='', layout='fruchtermanreingold', layout.weights=TRUE, plot.weights=TRUE,
			  node.size.title='', lab.use=NULL, ret.graph=FALSE){
    
    library(ggplot2)
    library(tidyverse)
    library(cowplot)
    library(ggnetwork)
    source('~/code/single_cell/colors.r')

    # plot network from graph or edgelist
    # -----------------------------------
    # graph = [m x n] matrix
    # edges = data.frame(1=source, 2=target, 3=weight, 4=color)
    # node_sizes = list(node_label = size)
    # node_colors = list(node_label = color)
    # edge_colors = [m x n] matrix
    # scale_nodes = c(min, max)
    # scale_edges = c(min, max)
    
    # convert graph to edgelist
    if(!is.null(graph)){
        graph = as.data.frame(as.matrix.data.frame(as.data.frame(graph)))
	edges = gather(graph %>% rownames_to_column('source'), target, weight, -source)
    }
    if(ncol(edges) == 3){colnames(edges) = c('source', 'target', 'weight')}
        
    # convert colors to edgelist
    if(!is.null(edge_colors)){
        edge_colors = as.data.frame(as.matrix.data.frame(as.data.frame(edge_colors)))
	edge_colors = gather(edge_colors %>% rownames_to_column('source'), target, color, -source)
	edges = merge(edges, edge_colors, by=c('source', 'target'))
    }
    if(! 'color' %in% colnames(edges)){edges$color = edge_color}
    
    # convert node_sizes to list
    nodes = sort(unique(c(as.character(edges$source), as.character(edges$target))))
    if(!is.null(node_sizes) & is.null(names(node_sizes))){names(node_sizes) = rownames(graph)} # assume graph
    if(is.null(node_sizes)){node_sizes = setNames(rep(1, length(nodes)), nodes)}
    
    # convert node_colors to list
    if(is.null(node_colors)){
        node_colors = setNames(rep(1, length(nodes)), nodes)
    } else {
        if(length(node_colors) == 1){
	    node_colors = setNames(rep(node_colors, length(nodes)), nodes)
	}
	if(is.null(names(node_colors))){names(node_colors) = rownames(graph)}
    }
    node_colors = as.factor(node_colors)
        
    # remove duplicate edges
    if(symm == TRUE){
        i = apply(edges[,1:2], 1, function(a) paste(sort(a), collapse=' '))
	edges = edges[!duplicated(i),]
    }
    
    # adjust edge weights
    edges = edges[edges$weight > 0,,drop=F]
    edges = edges[edges$source != edges$target,,drop=F]
    if(!is.null(scale_edges)){
        edges$weight = rescale_vector(edges$weight, target=scale_edges, abs=TRUE)
    }
    
    # convert edges to igraph object
    g = graph.data.frame(edges)
    
    # vertex attributes / igraph craziness
    nodes = V(g)$name
    
    # node sizes
    node_sizes = node_sizes[nodes]
    V(g)$node_size = node_sizes
    if(is.null(smin)){smin = min(node_sizes)}
    if(is.null(smax)){smax = max(node_sizes)}
    if(smin == smax){size_guide=FALSE} else {size_guide=TRUE}
    
    # node colors
    i = levels(node_colors)
    node_colors = node_colors[nodes]
    levels(node_colors) = i
    V(g)$node_color = as.character(node_colors)
    
    # igraph -> data.frame for ggplot2
    if(layout.weights == TRUE){
        net = ggnetwork(g, layout=layout, weights='weight', arrow.gap=1)
    } else {
        net = ggnetwork(g, layout=layout, arrow.gap=1)
    }
    
    if(nlevels(node_colors) > 0){
        net$node_color = factor(net$node_color, levels=levels(node_colors))
    }
        
    if(plot.weights == FALSE){
        net$weight[!is.na(net$weight)] = 1
    }
    return(net)
    # plot with ggnetwork
    if(is.null(out)){alpha = 1}
    
    if(length(unique(net$node_color)) > 1){
        p = ggplot(net, aes(x=x, y=y, xend=xend, yend=yend)) +
            geom_edges(color='#cccccc', size=na.omit(net$weight), alpha=alpha, curvature=curvature) +
            geom_nodes(aes(color=node_color, size=node_size)) +	    
            scale_color_manual(legend_title, breaks=levels(net$node_color), labels=levels(net$node_color), values=set.colors, drop=FALSE) +
	    scale_size(node.size.title, range=scale_nodes, limits=c(smin, smax), guide=size_guide)
    } else {
        p = ggplot(net, aes(x=x, y=y, xend=xend, yend=yend)) +
            geom_edges(color=na.omit(net$color), size=na.omit(net$weight), alpha=alpha, curvature=curvature) +
            geom_nodes(color='#0191C8', aes(size=node_size)) +
	    scale_size(node.size.title, range=scale_nodes, limits=c(smin, smax), guide=size_guide) +
	    scale_color_manual(guide=F)
    }
    
    if(is.null(lab.use)){
        lab.use = net$vertex.names
    }
    p = p + geom_nodetext_repel(aes(label=ifelse(vertex.names %in% lab.use, as.character(vertex.names), ''))) + theme_blank()
        
    if(ggtitle != ''){p = p + ggtitle(ggtitle)}
    
    if(do.legend == FALSE){p = p + theme(legend.position='none')}

    if(!is.null(out)){save_plot(p, file=out, nrow=nrow, ncol=ncol)}

    if(ret.graph == FALSE){p} else {list(p=p, g=g)}
}


nice_network = function(graph=NULL, edges=NULL, symm=TRUE, self=FALSE,
	                layout='fruchtermanreingold', layout_weights=TRUE,
			size_legend='nodes', size_trans='identity', # 'log2', etc
                        node_colors=NULL, node_sizes=NULL, node_scale=c(5,5), node_limits=NULL, node_alpha=1, node_pal=NULL, node_padding=.5, node_strokes=0.25,
			edge_colors=NULL, edge_alpha=1, edge_scale=c(.5,.5), edge_limits=NULL, edge_pal=NULL, color_sign=FALSE, edge_padding=.5, linewidth_pal=1:6,
			node_guide='legend', edge_guide='legend', size_guide='legend', linetype_guide='legend', do.legend=TRUE,
			node_title='', edge_title='', size_title='', linetype_title='',
			lab.use=NULL, node_lab.size=4, edge_lab.size=3, segment.color='grey',
			ggtitle='', drop_levels=FALSE, out=NULL, nrow=1, ncol=1, ret.data=FALSE,
			do.arrow=FALSE, arrow.gap=.015
			) {
    
    library(ggplot2)
    library(cowplot)    
    library(ggnetwork)
    library(igraph)
    library(ggrepel)
    library(tidyverse)
    source('~/code/single_cell/colors.r')
    
    # plot network from graph or edgelist
    # -----------------------------------
    # graph = [m x n] matrix
    # edges = data.frame(1=source, 2=target, 3=weight, 4=color, 5=size, 6=label, 7=linetype, ...)
    # lab.use = list of node/edge labels to use
    # node_sizes = list(label=size)
    # node_colors = list(label=color)
    # edge_colors = [m x n] matrix
    # node_scale = c(min, max)  * scales specify the minimum and maximium size of the plotting symbol
    # edge_scale = c(min, max) 
    # node_limits = c(min, max) * limits specify the minimum and maximum size of the value being plotted
    # edge_limits = c(min, max)
    
    # set arrow parameters
    # --------------------
    if(do.arrow == FALSE){
        arrow = NULL
	arrow.gap = 0
    } else {
        arrow = arrow(length=unit(6, "pt"), type="closed")
	arrow.gap = arrow.gap
    }
    
    # construct edgelist
    # ------------------
    if(!is.null(graph)){
        graph = as.data.frame(as.matrix.data.frame(as.data.frame(graph)))
	edges = gather(graph %>% rownames_to_column('source'), target, weight, -source)
    }
    if(! all(c('source', 'target', 'weight') %in% colnames(edges))){
        stop('check format')
    }
    
    # fix edgelist
    # ------------
    edges = as.data.table(edges)
    edges$source = as.character(edges$source)
    edges$target = as.character(edges$target)
    
    if(symm == TRUE){
        i = apply(edges[,.(source, target)], 1, function(a) paste(a, collapse=' '))
	edges = edges[!duplicated(i)]
    }
    if(self == FALSE){
        edges = edges[source != target]
    }
    edges = edges[abs(weight) > 0]
    
    # node attributes
    # ---------------
    node_labels = sort(unique(c(edges$source, edges$target)))

    node_sizes = unlist(node_sizes) # recently added
    if(is.null(node_sizes)){
        node_sizes = setNames(rep(1, length(node_labels)), node_labels)
    }
    if(is.null(names(node_sizes))){
        names(node_sizes) = node_labels
    }
    
    node_colors = unlist(node_colors) # recently added
    if(is.null(node_colors)){
        node_colors = setNames(rep(1, length(node_labels)), node_labels)
    }
    if(is.null(names(node_colors))){
        names(node_colors) = node_labels
    }
    if(is.null(levels(node_colors))){
        node_colors = as.factor(node_colors)
    }
    
    node_strokes = unlist(node_strokes) # recently added
    if(length(node_strokes) == 1){
        node_strokes = setNames(rep(node_strokes, length(node_labels)), node_labels)
    }
    
    
    # edge attributes
    # ---------------
    if(!is.null(edge_colors)){
        edge_colors = as.data.frame(as.matrix.data.frame(as.data.frame(edge_colors)))
	edge_colors = gather(edge_colors %>% rownames_to_column('source'), target, color, -source)
	edges = merge(edges, edge_colors, by=c('source', 'target'))
    }
    if(! 'color' %in% colnames(edges)){
        edges$color = 1
    }
    if(! 'label' %in% colnames(edges)){
        edges$label = ''
    }
    if(! 'linetype' %in% colnames(edges)){
        edges$linetype = 1
    }
    edges$color = as.factor(edges$color)
    edges$linetype = as.factor(edges$linetype)
    
    # fix edge weights
    # ----------------
    if(any(edges$weight < 0)){
        print('Warning: negative edges')
    }
    if(color_sign == TRUE){
        edges$color = factor(ifelse(edges$weight > 0, 'Positive', 'Negative'), levels=c('Positive', 'Negative'))
    }
    edges$weight = abs(edges$weight)
        
    # graph layout
    # ------------
    g = graph.data.frame(edges)
    node_order = V(g)$name
    V(g)$node_size = node_sizes[node_order]
    V(g)$node_color = as.character(node_colors[node_order])
    V(g)$node_stroke = as.numeric(node_strokes[node_order])
    layout_weights = ifelse(layout_weights == TRUE, 'weight', NULL)
    d = ggnetwork(g, layout=layout, weights=layout_weights, arrow.gap=arrow.gap)
    
    # scale sizes (after layout)
    # --------------------------
    i = (d$x == d$xend & d$y == d$yend)
    #d[!i, 'weight'] = rescale_vector(d[!i, 'weight'], edge_scale)
    if(is.null(node_limits)){
        node_limits = c(min(V(g)$node_size), max(V(g)$node_size))
    }
    
    # plot parameters
    # ---------------
    if(!is.null(lab.use)){
        d$vertex.names = ifelse(d$vertex.names %in% lab.use, as.character(d$vertex.names), '')
	d$label = ifelse(d$label %in% lab.use, as.character(d$label), '')
    }
    d$node_color = factor(d$node_color, levels=levels(node_colors))
    d$color = factor(d$color, levels=levels(edges$color))
    d$linetype = factor(d$linetype, levels=levels(edges$linetype))
    if(is.null(node_pal)){
        if(nlevels(d$node_color) == 1){node_pal = 'black'} else {node_pal = set.colors}
    }
    if(is.null(edge_pal)){
        if(nlevels(d$color) == 1){edge_pal = 'grey'} else {edge_pal = set.colors}
    }

    # self edges
    # ----------
    if(self == TRUE){
        self_nodes = edges[edges$source == edges$target, 'source']
	d$node_stroke = ifelse(i & d$source %in% self_nodes, d$node_stroke, 0)
    }
    
    # plot network
    # ------------
    p = ggplot(d)
    print(d)
    # scale nodes or edges (can't do both)
    if(size_legend == 'nodes'){
        p = p + geom_segment(data=d[!i,], aes(x=x, y=y, xend=xend, yend=yend, colour=color, linetype=linetype), size=d[!i,'weight'], alpha=edge_alpha, arrow=arrow) +
	    geom_point(data=d[i,], pch=21, aes(x=x, y=y, fill=node_color, size=node_size, stroke=node_stroke), alpha=node_alpha) +
	    scale_radius(size_title, limits=node_limits, range=node_scale, guide=size_guide, trans=size_trans)
    } else {
        p = p + geom_segment(data=d[!i,], aes(x=x, y=y, xend=xend, yend=yend, colour=color, linetype=linetype, size=weight), alpha=edge_alpha, arrow=arrow) +
	    geom_point(data=d[i,], pch=21, aes(x=x, y=y, fill=node_color, stroke=node_stroke), size=max(node_scale), alpha=node_alpha) +
	    scale_size(size_title, limits=edge_limits, range=edge_scale, guide=size_guide, trans=size_trans)
    }

    # construct the rest of the plot
    p = p + 
	geom_text_repel(data=d[i,], aes(x=x, y=y, label=vertex.names), box.padding=unit(node_padding, 'lines'), segment.color=segment.color, size=node_lab.size) +
	geom_text_repel(data=d[!i,], aes(x=(x+xend)/2, y=(y+yend)/2, label=label), box.padding=unit(edge_padding, 'lines'), segment.color=segment.color, size=edge_lab.size) +
	scale_fill_manual(node_title, breaks=levels(d$node_color), labels=levels(d$node_color), values=node_pal, guide=node_guide, drop=FALSE) +
	scale_colour_manual(edge_title, breaks=levels(d$color), labels=levels(d$color), values=edge_pal, guide=edge_guide, drop=FALSE) +
	scale_linetype_manual(linetype_title, guide=linetype_guide, drop=FALSE, values=linewidth_pal) +
	theme_blank() + ggtitle(ggtitle)
    
    # fix legends
    # -----------
    p = p + guides(
                fill = guide_legend(override.aes = list(size=.75*max(node_scale))),
		colour = guide_legend(override.aes = list(size=2))
		)
    if(do.legend == FALSE){
        p = p + theme(legend.position='none')
    }
    if(nlevels(d[i, 'node_color']) == 1){
        p = p + guides(fill=FALSE)
    }
    if(nlevels(d[!i, 'color']) == 1){
        p = p + guides(colour=FALSE)
    }
    if(size_legend == 'nodes'){
        if(min(d[i, 'node_size'], na.rm=T) == max(d[i, 'node_size'], na.rm=T)){p = p + guides(size=FALSE)}
    }
    if(nlevels(d[!i, 'linetype']) == 1){
        p = p + guides(linetype=FALSE)
    }
    p
    
    # write output
    # ------------
    if(!is.null(out)){
        save_plot(p, file=out, nrow=nrow, ncol=ncol)
    }
    if(ret.data == FALSE){
        p        
    } else {
        list(p=p, d=d)
    }
}
