// This is the script that contains all of the interactions for the main state estimation view.
// have a user-side cache to avoid having to re-make so many json calls
var cache = {};
cache.barplots = {};
cache.scatterplots = {};
cache.enrichr = {};
cache.cellmarker = {};

var currently_selected_cluster = 0;

var current_scatterplot_data = {};
var current_cluster_values = {};

var current_barplot_data = {};

var all_selected_clusters = [];
var current_selected_cells = [];
var currently_merging = false;

// criterion_template is used in custom_selections.js
var criterion_template = '';

// if this is true, then stop any update_barplot jobs.
var update_barplot_is_running = false;

// use jquery for ajax calls

function bind_click() {
    var plot = $('#means-scatter-plot')[0];
    plot.on('plotly_click', function(data) {
        var cluster = data.points[0].curveNumber;
        current_selected_cells = [];
        for (var i = 0; i<data.points.length; i++) {
            current_selected_cells.push(data.points[i].text);
        }
        // update selected cluster
        currently_selected_cluster = cluster;
        var top_or_bulk = $("#top-or-bulk").val();
        // set cluster 1 to be the selected value
        if (top_or_bulk.includes('pairwise')) {
            $('#barplot_cluster_select_1').val(String(data.points[0].curveNumber));
        }
        $('#cell_search_cluster').val(data.points[0].curveNumber);
        console.log(cluster);
        update_barplot(cluster);
        // update barplot
        var l = current_scatterplot_data.data[cluster].x.length;
        all_selected_clusters = [cluster];
        selection_string = "cluster " + String(cluster) + 
            " (" + String(l) + " cells), ";
        // set update area...
        $("#update-area").append("\nselected clusters: " + selection_string);
    });
}

// function called whenever a scatterplot is clicked or
// when the dropdown is changed...
function update_barplot(cluster_number) {
    if (update_barplot_is_running) {
        console.log('update_barplot is already running');
        return false;
    }
    update_barplot_is_running = true;
    var top_or_bulk = $("#top-or-bulk").val();
    var input_value = cluster_number;
    var num_genes = $("#num-genes").val();
    var cell_color = $("#cell-color").val();
    console.log('update_barplot');
    var cluster1 = $("#barplot_cluster_select_1").val();
    var cluster2 = $("#barplot_cluster_select_2").val();
    var selected_gene = $("#barplot_gene_select").val();
    var data = {"top_or_bulk": top_or_bulk,
               "input_value": input_value,
               "all_selected_clusters": all_selected_clusters,
               "num_genes": num_genes,
               "cell_color": cell_color,
               "cluster1": cluster1,
               "cluster2": cluster2,
               "selected_gene": selected_gene};
    if (top_or_bulk == 'double_pairs_comparison') {
        data['cluster1'] = $('#double_pair_cluster_1').val();
        data['cluster2'] = $('#double_pair_cluster_2').val();
        data['cluster3'] = $('#double_pair_cluster_3').val();
        data['cluster4'] = $('#double_pair_cluster_4').val();
    }
    if (top_or_bulk == 'violin') {
        data['violin_use_log'] = $('#violin_use_log').is(':checked') ? 1 : 0;
        data['violin_all_clusters'] = $('#violin_all_clusters').is(':checked') ? 1 : 0;
    }
    console.log(data);
    var key = JSON.stringify(data);
    $("#update-area").empty();
    $("#update-area").append('<span id="update_barplot">Updating barplot <img src="/static/ajax-loader.gif"/></span>');
    if (cache.barplots.hasOwnProperty(key)) {
        update_barplot_is_running = false;
        data = cache.barplots[key];
        if (data.data[0].type != 'histogram') {
            var gene_names = data.data[0].y;
            $('#top-genes-view').val(gene_names.join('\n'));
        }
        current_barplot_data = data;
        Plotly.newPlot("top-genes", data.data, data.layout, config={showSendToCloud: true});
        if (top_or_bulk == 'volcano_pairwise' || top_or_bulk == 'double_pairs_comparison') {
            var view = $('#top-genes')[0];
            view.on('plotly_selected', on_gene_select);
        }
        $('#update_barplot').remove();
        $("#update-area").append('Barplot updated');
        return true;
    }
    $.ajax({url: window.location.pathname + "/update_barplot",
        type: "POST",
        data: data,
    }).done(function(return_data) {
        update_barplot_is_running = false;
        if (return_data.startsWith('Error')) {
            $('#update_barplot').remove();
            $("#update-area").append(return_data);
            return false;
        }
        return_data = JSON.parse(return_data);
        cache.barplots[key] = return_data;
        if (return_data.data[0].type == 'bar') {
            var gene_names = return_data.data[0].y;
            $('#top-genes-view').val(gene_names.join('\n'));
            if (gene_names.length > 20) {
                return_data.layout.yaxis = {};
                return_data.layout.yaxis.showticklabels = false;
            }
            return_data.data[0].x.reverse();
            return_data.data[0].y.reverse();
        }
        current_barplot_data = return_data;
        Plotly.newPlot("top-genes", return_data.data, return_data.layout, config={showSendToCloud: true});
        if (top_or_bulk == 'volcano_pairwise' || top_or_bulk == 'double_pairs_comparison') {
            var view = $('#top-genes')[0];
            view.on('plotly_selected', on_gene_select);
        }
        $('#update_barplot').remove();
        $("#update-area").append('Barplot updated');
    });
    return true;
}

// TODO: update the barplot for gene selection
function update_barplot_genes() {
}

// saves the scatterplot as an svg
function download_scatterplot_svg() {
    var plot = JSON.parse(JSON.stringify(current_scatterplot_data));
    for (var i = 0; i < plot.data.length; i++) {
        if (plot.data[i].type == 'scattergl') {
            plot.data[i].type = 'scatter';
        }
    }
    Plotly.downloadImage(plot, {format: 'svg'});
}

// saves the barplot as an svg
function download_barplot_svg() {
    var plot = $('#top-genes')[0];
    Plotly.downloadImage(plot, {format: 'svg'});
}

// TODO: download scatterplot/barplot data as a tab-separated table
function download_scatterplot_data() {
    // TODO: look at the different data types, convert to tsv
    var plot_type = $('input[name="scatter-type"]:checked').val();
    var data = current_scatterplot_data.data;
    // scatterplot data
    if (plot_type == 'Baseline' || plot_type == 'Cells' || plot_type == 'Means' || plot_type == 'Genes') {
        var output_string = 'cluster_name\tid\tx\ty\n';
        for (var track_index in data) {
            var track = data[track_index];
            var track_length = track.x.length;
            for (var i = 0; i < track_length; i++) {
                output_string += track.name + '\t' + track.text[i] + '\t' + track.x[i] + '\t' + track.y[i] + '\n';
            }
        }
        var output = "data:application/octet-stream," + encodeURI(output_string);
        location.href = output;
    } else if (plot_type == 'Dendrogram') {
        var heatmap_data = data[data.length - 1];
        var cluster_names = current_scatterplot_data.layout.xaxis.ticktext;
        var gene_names = current_scatterplot_data.layout.yaxis.ticktext;
        var output_string = 'Gene_name\t';
        for (var i in cluster_names) { 
            output_string += cluster_names[i] + '\t';
        }
        output_string += '\n';
        for (var i in gene_names) {
            var color_data = heatmap_data.z[i];
            output_string += gene_names[i] + '\t';
            for (var j in cluster_names) {
                output_string += color_data[j] + '\t';
            }
            output_string += '\n'
        }
        var output = "data:application/octet-stream," + encodeURI(output_string);
        location.href = output;
    } else if (plot_type == 'Cluster_heatmap' || plot_type == 'Gene_heatmap' || plot_type == 'Diffcorr_heatmap' || plot_type == 'Correlation_heatmap') {
        var x_names = data[0].x;
        var y_names = data[0].y;
        var output_string = 'y\t';
        for (var i in x_names) {
            output_string += x_names[i] + '\t';
        }
        output_string += '\n';
        for (var i in y_names) {
            var color_data = data[0].z[i];
            output_string += y_names[i] + '\t';
            for (var j in x_names) {
                output_string += color_data[j] + '\t';
            }
            output_string += '\n'
        }
        var output = "data:application/octet-stream," + encodeURI(output_string);
        location.href = output;
    }
}

function download_barplot_data() {
    var data = current_barplot_data.data;
}

// toggle special areas for the scatterplot
function toggle_scatterplot_type() {
    var plot_type = $('input[name="scatter-type"]:checked').val();
    if (plot_type == 'Cluster_heatmap') {
        $('#heatmap-names-area').toggle(true);
        $('#dendrogram-names-area').toggle(false);
        $('#gene-heatmap-options-area').toggle(false);
        $('#diffcorr-heatmap-options-area').toggle(false);
    } else if (plot_type == 'Dendrogram') {
        $('#heatmap-names-area').toggle(false);
        $('#dendrogram-names-area').toggle(true);
        $('#gene-heatmap-options-area').toggle(false);
        $('#diffcorr-heatmap-options-area').toggle(false);
    } else if(plot_type == 'Gene_heatmap') {
        $('#heatmap-names-area').toggle(false);
        $('#dendrogram-names-area').toggle(false);
        $('#gene-heatmap-options-area').toggle(true);
        $('#diffcorr-heatmap-options-area').toggle(false);
    } else if(plot_type == 'Diffcorr_heatmap') {
        $('#heatmap-names-area').toggle(false);
        $('#dendrogram-names-area').toggle(false);
        $('#gene-heatmap-options-area').toggle(false);
        $('#diffcorr-heatmap-options-area').toggle(true);
    } else {
        $('#heatmap-names-area').toggle(false);
        $('#dendrogram-names-area').toggle(false);
        $('#gene-heatmap-options-area').toggle(false);
        $('#diffcorr-heatmap-options-area').toggle(false);
    }
}

// this function is called on startup, and whenever the radio buttons
// corresponding to different input types are clicked.
function update_scatterplot() {
    var plot_type = $('input[name="scatter-type"]:checked').val();
    var cell_color = $("#cell-color").val();
    var gene_name = $('#gene_name').val();
    console.log(plot_type);
    console.log(cell_color);
    var use_mw = false;
    if (cell_color == "gene") {
        use_mw = $('#use_mw_gene').val();
    }
    var cluster = $('#cluster_input').val();
    var upload_data = {"scatter_type": plot_type, "cell_color": cell_color,
               "gene_name": gene_name,
               "use_mw": use_mw,
               "cluster_input": cluster,
    };
    if (plot_type == 'Cluster_heatmap') {
        upload_data['heatmap_cluster_name_1'] = $('#heatmap_cluster_name_1').val();
        upload_data['heatmap_cluster_name_2'] = $('#heatmap_cluster_name_2').val();
    }
    if (plot_type == 'Dendrogram') {
        var dendrogram_genes = $('#dendrogram_genes').val();
        upload_data['dendrogram_genes'] = dendrogram_genes;
        upload_data['dendrogram_use_log'] = $('#dendrogram_use_log').is(':checked') ? 1 : 0;
        upload_data['dendrogram_normalize'] = $('#dendrogram_normalize').is(':checked') ? 1 : 0;
    }
    if (plot_type == 'Gene_heatmap') {
        upload_data['heatmap_genes_1'] = $('#heatmap_genes_1').val();
        upload_data['heatmap_genes_2'] = $('#heatmap_genes_2').val();
        upload_data['gene_heatmap_cluster'] = $('#gene_heatmap_cluster').val();
    }
    if (plot_type == 'Diffcorr_heatmap') {
        upload_data['diffcorr_genes_1'] = $('#diffcorr_genes_1').val();
        upload_data['diffcorr_genes_2'] = $('#diffcorr_genes_2').val();
        upload_data['diffcorr_cluster_1'] = $('#diffcorr_cluster_1').val();
        upload_data['diffcorr_cluster_2'] = $('#diffcorr_cluster_2').val();
        upload_data['diffcorr_value'] = $('#diffcorr_value').val();
    }
    $("#update-area").empty();
    $("#update-area").append('Updating scatterplot <img src="/static/ajax-loader.gif"/>');
    var key = JSON.stringify(upload_data);
    // if the plot parameters have been used before, we retrieve them from the cache...
    if (cache.scatterplots.hasOwnProperty(key)) {
        var data = cache.scatterplots[key];
        Plotly.newPlot("means-scatter-plot", data.data, data.layout, config={showSendToCloud:true});
        current_scatterplot_data = data;
        if (data.data.length > 1) {
            update_cluster_selections();
        }
        bind_click();
        bind_select();
        $("#update-area").empty();
        $("#update-area").append('Scatterplot updated');
        return true;
    }
    $.ajax({url: window.location.pathname + "/update_scatterplot",
        type: "POST",
        data: upload_data,
    }).done(function(data) {
        if (data.startsWith('Error')) {
            $("#update-area").empty();
            $("#update-area").append(data);
            return false;
        }
        data = JSON.parse(data);
        cache.scatterplots[key] = data;
        Plotly.newPlot("means-scatter-plot", data.data, data.layout, config={showSendToCloud:true});
        current_scatterplot_data = data;
        bind_click();
        bind_select();
        $("#update-area").empty();
        $("#update-area").append('Scatterplot updated');
        // update selections that require cluster names
        if (data.data.length > 1) {
            update_cluster_selections();
        }
    });
    return true;
}

// updates all fields that involve selecting clusters.
// This is called whenever the scatterplot changes...
function update_cluster_selections() {
    var s1 = $('#barplot_cluster_select_1');
    var s2 = $('#barplot_cluster_select_2');
    var s3 = $('#gene_heatmap_cluster');
    var s4 = $('#diffcorr_cluster_1');
    var s5 = $('#diffcorr_cluster_2');
    var s6 = $('#cell_search_cluster');
    var cluster_selects = [s1, s2, s3, s4, s5, s6, $('#double_pair_cluster_1'), $('#double_pair_cluster_2'), $('#double_pair_cluster_3'), $('#double_pair_cluster_4')];
    for (var s in cluster_selects) {
        var select = cluster_selects[s];
        select.empty();
        if (s == 2) {
            select.append($('<option>').attr('value', 'all').text('All cells'));
        }
        for (var i in current_scatterplot_data.data) {
            var value = current_scatterplot_data.data[i];
            select.append($("<option>").attr('value', i).text(value.name));
        }
    }
}

function set_enrichr_results(data, query) {
    // create table
    // TODO: make it collapsible
    var results_view = $('#' + query + '_results');
    results_view.empty();
    var table = $('<table>');
    var header = $('<tr>');
    for (var j = 0; j<data[0].length; j++) {
        header.append('<th>'+data[0][j]+'</th>');
    }
    table.append(header);
    for (var i = 1; i<data.length; i++) {
        var row = $('<tr>');
        for (var j = 0; j<data[i].length; j++) {
            row.append('<td>'+data[i][j]+'</td>');
        }
        table.append(row);
    }
    results_view.append(table);
}

function update_gene_query(query) {
    // query is a string that can be 'enrichr', 'cellmarker', etc...
    var input_array = $('#gene-set-query-form').serializeArray();
    var data = {};
    $(input_array).each(function(index, obj){
        data[obj.name] = obj.value;
    });
    console.log('update_gene_query');
    console.log(data);
    var key = JSON.stringify(data);
    if (cache.enrichr.hasOwnProperty(key)) {
        set_enrichr_results(cache.enrichr[key], query);
        return true;
    }
    var results_view = $('#' + query + '_results');
    var update_url = '/update_' + query;
    results_view.empty();
    results_view.append("<br>" + "Query in progress..." + '<img src="/static/ajax-loader.gif"/>');
    $.ajax({url: window.location.pathname + update_url,
        type: "POST",
        data: data,
    }).done(function(results) {
        results_view.empty();
        if (results.startsWith('Error')) {
            results_view.append(results);
        } else {
            results = JSON.parse(results);
            cache.enrichr[key] = results;
            set_enrichr_results(results, query);
        }
    });
    return true;
}

function bind_select() {
   var plot = $('#means-scatter-plot')[0];
   plot.on('plotly_selected', on_select);
}

// handler for scatterplot selection events
function on_select(data) {
    console.log('on_select');
    console.log(data);
    if (!data) {
        return false;
    }
    var cell_color = $("#cell-color").val();
    var clusters = new Set([]);
    // identify selected cells
    current_selected_cells = [];
    for (var i = 0; i<data.points.length; i++) {
        clusters.add(data.points[i].curveNumber);
        current_selected_cells.push(data.points[i].text);
    }
    // TODO: this uses uncurl clusters, not the currently selected cell labels
    // find cluster lengths
    var cluster_lengths = [];
    all_selected_clusters = [];
    var selection_string = "";
    clusters.forEach(function(x) {
        var l = current_scatterplot_data.data[x].x.length;
        cluster_lengths.push(l);
        all_selected_clusters.push(x);
        selection_string += "cluster " + String(x) + 
            " (" + String(l) + " cells), ";
    });
    // set update area...
    $("#update-area").empty();
    $("#update-area").append("selected clusters: " + selection_string + "<br>");
    $("#update-area").append("Number of selected cells: " + current_selected_cells.length);
    // update custom criteria cell selection
    set_selection_target_to_selected_cells();
}

// function for handling selections on a gene scatter plot
// send it to the bottom view
function on_gene_select(data) {
    console.log('on_gene_select');
    if (!data) {
        return false;
    }
    var gene_names = [];
    for (var i = 0; i < data.points.length; i++) {
        gene_names.push(data.points[i].text);
    }
    // push gene names to the box
    $('#top-genes-view').val(gene_names.join('\n'));
}

// function for handling clicks on the gene scatter plot
// select all genes within that cluster and send them to the query box.
function on_gene_click(data) {
    // TODO
    console.log('on_gene_click');
    if (!data) {
        return false;
    }
}

// Gets basic cell info - read counts and gene counts
function get_cell_info() {
    var selected_clusters = all_selected_clusters;
    var selected_cells = current_selected_cells;
    // can't send an array... have to convert to string
    // TODO: get color map
    var upload_data = {
        selected_clusters: String(selected_clusters),
        selected_cells: String(selected_cells),
        color_map: $('#cell-color').val()
    };
    $("#update-area").empty();
    $("#update-area").append("<br>" + "Query in progress..." + '<img src="/static/ajax-loader.gif"/>');
    console.log(upload_data);
    $.ajax({url: window.location.pathname + "/cell_info",
        data: upload_data,
        method: 'POST',
    }).done(function(data) {
        var results = JSON.parse(data);
        var selection_string = "";
        selected_clusters.forEach(function(x) {
            var l = current_scatterplot_data.data[x].x.length;
            selection_string += "cluster " + String(x) + 
                " (" + String(l) + " cells), ";
        });
        var cell_string = "";
        for (var i = 0; i < selected_cells.length; i++) {
            cell_string += "cell " + String(selected_cells[i]) + " (" + String(results.gene_counts[i])
                +  " genes, " + String(results.read_counts[i]) + " reads), ";
        }
        $("#update-area").empty();
        $("#update-area").append("Selected cells: " + cell_string + "<br>");
        $("#update-area").append("Selected clusters: " + selection_string + "<br>");
        $("#update-area").append("Median read count for cluster: " + results.cluster_reads + "<br>");
        $("#update-area").append("Median gene count for cluster: " + results.cluster_genes + "<br>");
        $("#update-area").append("Total nonzero gene count for cluster: " + results.total_gene_count + "<br>");
    });
}

function split_or_merge_cluster(split_or_merge, cells_or_clusters) {
    if (currently_merging) {
        window.alert('Already running ' + split_or_merge + ' cluster');
        return false;
    }
    var result = window.confirm("Warning: splitting or merging might take a while, depending on the size of the dataset. Are you sure you wish to run this operation?");
    if (result == false) {
        return false;
    }
    $('.overlay').show();
    var selected_clusters = all_selected_clusters;
    if (cells_or_clusters == "cells") {
        selected_clusters = current_selected_cells;
    }
    if (all_selected_clusters.length == 0) {
        window.alert('Warning: no selected clusters.');
        $('.overlay').hide();
        return false;
    }
    if (split_or_merge == "merge" && all_selected_clusters.length == 1) {
        window.alert('Warning: cannot merge a single cluster.');
        $('.overlay').hide();
        return false;
    }
    // lengths of each of the selected clusters
    console.log(split_or_merge);
    console.log(selected_clusters);
    // make the whole page block or something like that
    $("#update-area").append("<br>" + split_or_merge + " clusters in progress... (re-running UNCURL, recalculating differentially expressed genes) " + '<img src="/static/ajax-loader.gif"/>');
    currently_merging = true;
    // make some indication that split/merge has been called.
    $.ajax({url: window.location.pathname + "/split_or_merge_cluster",
        data: {'split_or_merge': split_or_merge,
               'selected_clusters': selected_clusters.join(',')},
        method: 'POST',
    }).done(function(data) {
        currently_merging = false;
        $('.overlay').hide();
        if (data.startsWith('Error')) {
            $('#update-area').empty();
            $('#update-area').append(data);
        } else {
            $('#update-area').empty();
            $('#update-area').append(data);
            // reload page?
            cache.barplots = {};
            cache.scatterplots = {};
            update_scatterplot();
            update_barplot(0);
            //location.reload(true);
            get_history();
        }
    });
}

// re-run pipeline
function delete_rerun() {
    var result = window.confirm("Do you wish to delete the current results for this dataset and re-run the entire pipeline?");
    if (result == false) {
        return false;
    }
    window.location.href = window.location.pathname + "/delete_rerun";
}


// re-run pipeline
function run_batch_correction() {
    var cell_color = $("#cell-color").val();
    var result = window.confirm("Do you wish to run batch effect correction on the colormap '" + cell_color + "'? This will re-run the entire pipeline, deleting the current dataset.");
    if (result == false) {
        return false;
    }
    $('.overlay').show();
    $.ajax({url: window.location.pathname + "/run_batch_correction",
        method: 'POST',
        data: {'colormap': cell_color,},
    }).done(function(data) {
        $('.overlay').hide();
        if (data.startsWith('Error')) {
            $('#update-area').empty();
            $('#update-area').append(data);
        } else {
            $('#update-area').empty();
            $('#update-area').append(data);
            // reload page?
            cache.barplots = {};
            cache.scatterplots = {};
            update_scatterplot();
            update_barplot(0);
            //location.reload(true);
            get_history();
        }
    });
}

// copies data to a new user_id
function copy_data() {
    var result = window.confirm("Do you wish to copy this dataset to a new user id?");
    if (result == false) {
        return false;
    }
    $('#update-area').empty();
    $("#update-area").append("Copying data and results " + '<img src="/static/ajax-loader.gif"/>');
    $.ajax({url: window.location.pathname + "/copy_dataset",
        method: 'POST'
    }).done(function(data) {
        if (data.startsWith('Error')) {
            $('#update-area').empty();
            $('#update-area').append(data);
        } else {
            var new_user_id = data;
            $('#update-area').empty();
            $('#update-area').append('New results page: <a href="/user/' + new_user_id + '/view">' + new_user_id + '</a>');
        }
    });
}


// Starts a new uncurl sub-analysis...
function subset_cells(is_cells) {
    var result = window.confirm("Do you wish to run a new analysis on the selected cells?");
    if (result == false) {
        return false;
    }
    $('#update-area').empty();
    $("#update-area").append("Copying data, generating summary " + '<img src="/static/ajax-loader.gif"/>');
    var cell_ids = current_selected_cells;
    if (!is_cells) {
        cell_ids = all_selected_clusters;
    }
    console.log(is_cells);
    console.log(cell_ids);
    $.ajax({url: window.location.pathname + "/subset",
        data: {
            'is_cells': is_cells ? 1 : 0,
            'cell_ids': cell_ids.join(','),
            'color_map': $('#cell-color').val(),
        },
        method: 'POST'
    }).done(function(data) {
        if (data.startsWith('Error')) {
            $('#update-area').empty();
            $('#update-area').append(data);
        } else {
            var new_user_id = data;
            $('#update-area').empty();
            $('#update-area').append('New results page: <a href="/state_estimation/results/' + new_user_id + '">' + new_user_id + '</a>');
        }
    });
}

function toggle_reanalyze_area(value) {
    if (value == 'reanalyze') {
        $('#reanalyze-area').toggle(true);
        $('#recluster_area').toggle(false);
    } else if (value == 'recluster') {
        $('#reanalyze-area').toggle(false);
        $('#recluster_area').toggle(true);
    }
}

// toggles visibility of the bottom-left gene query view
// 'toggle' is a jquery method that changes the visibility of the element.
function toggle_query_visibility() {
    var views = ['cellmarker_view', 'enrichr_view', 'cellmesh_view', 'cellmesh_anatomy_view', 'go_view', 'subtiwiki_view', 'kegg_view'];
    var value = $('#database-select').val();
    var view_select = value + '_view';
    for (var i in views) {
        if (views[i] == view_select) {
            $('#'+views[i]).toggle(true);
        } else {
            $('#'+views[i]).toggle(false);
        }
    }
}

// re-runs the clustering process, perhaps with a different algorithm.
function rerun_clustering() {
    var result = window.confirm('Warning: this may take a while. Continue?');
    if (result == false) {
        return 0;
    }
    var clustering_method = $('#clustering_method_select').val();
    console.log(clustering_method);
    $("#update-area").empty();
    $('.overlay').show();
    $("#update-area").append('Re-clustering + recalculating differential expression... <img src="/static/ajax-loader.gif"/>');
    $.ajax({url: window.location.pathname + "/recluster",
        data: {
            'clustering_method': clustering_method
        },
        method: 'POST'
    }).done(function(data) {
        $('.overlay').hide();
        if (data.startsWith('Error')) {
            $('#update-area').empty();
            $('#update-area').append(data);
        } else {
            var new_user_id = data;
            $('#update-area').empty();
            $('#update-area').append('Finished re-clustering.');
            // clear cache
            cache.barplots = {};
            cache.scatterplots = {};
            update_scatterplot();
        }
    });
}

function submit_db_query() {
    var form_data = $('#cell_search_form').serializeArray();
    var data = {};
    $(form_data).each(function(index, obj){
        data[obj.name] = obj.value;
    });
    var cell_color = $('#cell-color').val();
    data['cell_color'] = cell_color;
    $("#update-area").empty();
    $('#update-area').append('Cell search query <img src="/static/ajax-loader.gif"/>');
    $.ajax({url: window.location.pathname + '/db_query',
        data: data,
        method: 'POST'
    }).done(function(data) {
        if (data.startsWith('Error')) {
            $("#update-area").empty();
            $("#update-area").append(data);
            return false;
        }
        // set results from cell search...
        $("#update-area").empty();
        $("#update-area").append('Completed cell search query');
        var results = JSON.parse(data);
        set_enrichr_results(results, 'cell_search');
    });
}

// called whenever cell color is changed.
function on_cell_color_change() {
    var cell_color = $("#cell-color").val();
    if (cell_color == "gene") {
        var gene_name = $('#gene_name').val();
        if (gene_name.length > 0) {
            update_scatterplot();
        }
        $('#gene-name-area').css('display', 'block');
        $('#color-track-upload-area').css('display', 'none');
        $('#cluster-select-area').css('display', 'none');
        $('#custom_color_map_area').css('display', 'none');
    } else if (cell_color == "new") {
        $('#color-track-upload-area').css('display', 'block');
        $('#gene-name-area').css('display', 'none');
        $('#cluster-select-area').css('display', 'none');
        $('#custom_color_map_area').css('display', 'none');
    } else if (cell_color == "weights") {
        $('#cluster-select-area').css('display', 'block');
        $('#gene-name-area').css('display', 'none');
        $('#color-track-upload-area').css('display', 'none');
        $('#custom_color_map_area').css('display', 'none');
    } else if (cell_color == 'custom') {
        $('#gene-name-area').css('display', 'none');
        $('#color-track-upload-area').css('display', 'none');
        $('#cluster-select-area').css('display', 'none');
        $('#custom_color_map_area').css('display', 'block');
        add_custom_colormap();
    } else {
        $('#gene-name-area').css('display', 'none');
        $('#color-track-upload-area').css('display', 'none');
        $('#cluster-select-area').css('display', 'none');
        update_scatterplot();
        // this function checks if the label scheme is a custom label scheme,
        // and sends a json query to get the custom label schemes.
        get_custom_colormap();
    }
}

function on_top_or_bulk_change() {
    var option = $('#top-or-bulk').val();
    // if the selected option has "pairwise":
    if (option.includes('pairwise')) {
        $('#barplot-cluster-select-view').css('display', 'block');
        $('#barplot_double_pair_cluster_select').css('display', 'none');
        $('#barplot_violin_plot_options').css('display', 'none');
    } else {
        $('#barplot-cluster-select-view').css('display', 'none');
        if (option == "violin") {
            $('#barplot_double_pair_cluster_select').css('display', 'none');
            $('#barplot_violin_plot_options').css('display', 'block');
        } else if (option == "double_pairs_comparison") {
            $('#barplot_double_pair_cluster_select').css('display', 'block');
            $('#barplot_violin_plot_options').css('display', 'none');
        } else {
            $('#barplot_double_pair_cluster_select').css('display', 'none');
            $('#barplot_violin_plot_options').css('display', 'none');
            update_barplot(currently_selected_cluster);
        }
    }
}

// this function creates a restore_history function for the action_id.
// Used as an event.
function create_restore_history_function(action_id) {
    return function() {
        var result = window.confirm("Warning: this operation might take a while, and changes could be lost. Are you sure you wish to run this operation?");
        if (result == false) {
            return false;
        }
        $('.overlay').show();
        $.ajax({url: window.location.pathname + '/restore_history/' + action_id,
            method: 'GET'
        }).done(function(data) {
            $('.overlay').hide();
            if (data.startsWith('Error')) {
                $("#update-area").empty();
                $("#update-area").append(data);
                return false;
            }
            // TODO: refresh the screen/reload the page?
            window.location.reload();
        });
    };
}

function get_history() {
    // fills the history pane
    console.log('get_history');
    $.ajax({url: window.location.pathname + '/history',
        method: 'GET'
    }).done(function(data) {
        if (data.startsWith('Error')) {
            $("#update-area").empty();
            $("#update-area").append(data);
            return false;
        }
        var results = JSON.parse(data);
        $("#history_pane").empty();
        var tab = $('<table>');
        var el = $('<tr><td>Action</td><td>Time</td><td>Restore</td></tr>');
        tab.append(el);
        // TODO: link for "restoring", don't need ids (entry 1)
        for (var i = 0; i < results.length; i++) {
            var x = results[i];
            el = $('<tr>');
            var action = $('<td>').text(x[0]);
            el.append(action);
            var date = $('<td>').text(x[2]);
            el.append(date);
            var is_restorable = x[3];
            var id = x[1];
            if (is_restorable) {
                var restore_link = $('<button>');
                restore_link = restore_link.text('Restore previous');
                restore_link.click(create_restore_history_function(id));
                var restore_entry = $('<td>');
                restore_entry.append(restore_link);
                el.append(restore_entry);
            }
            tab.append(el);
        }
        $('#history_pane').append(tab);
    });
}


window.onload = function() {
    // activate tooltips
    $('[data-toggle="tooltip"]').tooltip();
    // bind radio buttons to update_scatterplot
    $('input[name="scatter-type"]').change(function() {
        update_scatterplot();
    });
    $('#cell-color').change(on_cell_color_change);

    update_scatterplot();

    // bind dropdowns to update_barplot
    $('#top-or-bulk').change(on_top_or_bulk_change);
    $('#num-genes').change(function() {
        update_barplot(currently_selected_cluster);
    });
    update_barplot(currently_selected_cluster);
    on_top_or_bulk_change();

    // bind update enrichr button to update_enrichr
    $('#enrichr-submit').click(function() {
        update_gene_query('enrichr');
    });

    $('#delete_rerun').click(function() {
        delete_rerun();
    });

    $('#copy').click(function() {
        copy_data();
    });

    $('#subset_cells').click(function() {
        // subset cells
        subset_cells(true);
    });

    $('#subset_clusters').click(function() {
        // subset clusters
        subset_cells(false);
    });

    criterion_template = document.getElementById('custom_selection_criterion-1').outerHTML;
    toggle_query_visibility();
};

