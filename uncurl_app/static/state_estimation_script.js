// This is the script that contains all of the interactions for the main state estimation view.
// have a user-side cache to avoid having to re-make so many json calls
var cache = {};
cache.barplots = {};
cache.scatterplots = {};
cache.enrichr = {};
cache.cellmarker = {};

var currently_selected_cluster = 0;

var current_scatterplot_data = {};

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
        //$("#update-area").empty();
        $("#update-area").append("\nselected clusters: " + selection_string);
    });
};

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
    var key = top_or_bulk + String(input_value) + ' ' + String(num_genes) + ' ' + String(cluster1) + ' ' + String(cluster2) + ' ' + String(cell_color) + ' ' + selected_gene;
    $("#update-area").empty();
    $("#update-area").append('Updating barplot <img src="/static/ajax-loader.gif"/>');
    if (cache.barplots.hasOwnProperty(key)) {
        update_barplot_is_running = false;
        var data = cache.barplots[key];
        if (data.data[0].type != 'histogram') {
            var gene_names = data.data[0].y;
            $('#top-genes-view').val(gene_names.join('\n'));
        }
        Plotly.newPlot("top-genes", data.data, data.layout);
        if (top_or_bulk == 'volcano_pairwise') {
            var view = $('#top-genes')[0];
            view.on('plotly_selected', on_gene_select);
        }
        $("#update-area").empty();
        $("#update-area").append('Barplot updated');
        return true;
    }
    $.ajax({url: window.location.pathname + "/update_barplot",
        type: "POST",
        data: {"top_or_bulk": top_or_bulk,
               "input_value": input_value,
               "num_genes": num_genes,
               "cell_color": cell_color,
               "cluster1": cluster1,
               "cluster2": cluster2,
               "selected_gene": selected_gene},
    }).done(function(data) {
        update_barplot_is_running = false;
        if (data.startsWith('Error')) {
            $("#update-area").empty();
            $("#update-area").append(data);
            return false;
        }
        data = JSON.parse(data);
        cache.barplots[key] = data;
        if (data.data[0].type == 'bar') {
            var gene_names = data.data[0].y;
            $('#top-genes-view').val(gene_names.join('\n'));
            if (gene_names.length > 20) {
                data.layout.yaxis = {};
                data.layout.yaxis.showticklabels = false;
            }
            data.data[0].x.reverse();
            data.data[0].y.reverse();
        }
        Plotly.newPlot("top-genes", data.data, data.layout);
        if (top_or_bulk == 'volcano_pairwise') {
            var view = $('#top-genes')[0];
            view.on('plotly_selected', on_gene_select);
        }
        $("#update-area").empty();
        $("#update-area").append('Barplot updated');
    });
    return true;
};

// TODO: update the barplot for gene selection
function update_barplot_genes() {
};

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
        var select = $('#gene_heatmap_cluster');
        // clear select children
        select.empty();
        select.append($('<option>').attr('value', 'all').text('All cells'));
        // get all cluster values
        for (var i in current_scatterplot_data.data) {
            var value = current_scatterplot_data.data[i];
            console.log('cluster ' + i);
            select.append($("<option>").attr('value', i).text(value.name));
        }
    } else if(plot_type == 'Diffcorr_heatmap') {
        $('#heatmap-names-area').toggle(false);
        $('#dendrogram-names-area').toggle(false);
        $('#gene-heatmap-options-area').toggle(false);
        $('#diffcorr-heatmap-options-area').toggle(true);
        var select = $('#diffcorr_cluster_1');
        var select2 = $('#diffcorr_cluster_2');
        // clear select children
        select.empty();
        select2.empty();
        // get all cluster values
        for (var i in current_scatterplot_data.data) {
            var value = current_scatterplot_data.data[i];
            console.log('cluster ' + i);
            select.append($("<option>").attr('value', i).text(value.name));
            select2.append($("<option>").attr('value', i).text(value.name));
        }
    } 
    else {
        $('#heatmap-names-area').toggle(false);
        $('#dendrogram-names-area').toggle(false);
        $('#gene-heatmap-options-area').toggle(false);
        $('#diffcorr-heatmap-options-area').toggle(false);
    }
};

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
        Plotly.newPlot("means-scatter-plot", data.data, data.layout);
        current_scatterplot_data = data;
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
        Plotly.newPlot("means-scatter-plot", data.data, data.layout);
        current_scatterplot_data = data;
        bind_click();
        bind_select();
        $("#update-area").empty();
        $("#update-area").append('Scatterplot updated');
        // update selections that require cluster names
        var cluster_select = $('#cell_search_cluster');
        cluster_select.empty();
        for (var i in current_scatterplot_data.data) {
            var value = current_scatterplot_data.data[i];
            cluster_select.append($("<option>").attr('value', i).text(value.name));
        }
    });
    return true;
};

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
};

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
};

function bind_select() {
   var plot = $('#means-scatter-plot')[0];
   plot.on('plotly_selected', on_select);
};

// handler for scatterplot selection events
function on_select(data) {
    console.log('on_select');
    if (!data) {
        return false;
    }
    var clusters = new Set([]);
    // identify selected cells
    current_selected_cells = [];
    for (var i = 0; i<data.points.length; i++) {
        clusters.add(data.points[i].curveNumber);
        current_selected_cells.push(data.points[i].text);
    }
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
};

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
};

// function for handling clicks on the gene scatter plot
// select all genes within that cluster and send them to the query box.
function on_gene_click(data) {
    // TODO
    console.log('on_gene_click');
    if (!data) {
        return false;
    }
};

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
        };
        $("#update-area").empty();
        $("#update-area").append("Selected cells: " + cell_string + "<br>");
        $("#update-area").append("Selected clusters: " + selection_string + "<br>");
        $("#update-area").append("Median read count for cluster: " + results.cluster_reads + "<br>");
        $("#update-area").append("Median gene count for cluster: " + results.cluster_genes + "<br>");
    });
};

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
        return false;
    }
    if (split_or_merge == "merge" && all_selected_clusters.length == 1) {
        window.alert('Warning: cannot merge a single cluster.');
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
        }
    });
};

// re-run pipeline
function delete_rerun() {
    var result = window.confirm("Do you wish to delete the current results for this dataset and re-run the entire pipeline?");
    if (result == false) {
        return false;
    }
    window.location.href = window.location.pathname + "/delete_rerun";
};

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
};


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
};

function toggle_reanalyze_area(value) {
    if (value == 'reanalyze') {
        $('#reanalyze-area').toggle(true);
        $('#recluster_area').toggle(false);
    } else if (value == 'recluster') {
        $('#reanalyze-area').toggle(false);
        $('#recluster_area').toggle(true);
    }
};

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
};

// re-runs the clustering process, perhaps with a different algorithm.
function rerun_clustering() {
    var result = window.confirm('Warning: this may take a while. Continue?');
    if (result == false) {
        return 0;
    }
    var clustering_method = $('#clustering_method_select').val();
    console.log(clustering_method);
    $("#update-area").empty();
    $("#update-area").append('Re-clustering + recalculating differential expression... <img src="/static/ajax-loader.gif"/>');
    $.ajax({url: window.location.pathname + "/recluster",
        data: {
            'clustering_method': clustering_method
        },
        method: 'POST'
    }).done(function(data) {
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
};

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
};

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
};


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
    $('#top-or-bulk').change(function() {
        var option = $('#top-or-bulk').val();
        // if the selected option has "pairwise":
        if (option.includes('pairwise')) {
            $('#barplot-cluster-select-view').css('display', 'block');
            // set the select values for the current color scheme...
            var select_1 = $('#barplot_cluster_select_1');
            var select_2 = $('#barplot_cluster_select_2');
            // clear select children
            select_1.empty();
            select_2.empty();
            // get all cluster values
            for (var i in current_scatterplot_data.data) {
                var value = current_scatterplot_data.data[i];
                console.log('cluster ' + i);
                select_1.append($("<option>").attr('value', i).text(value.name));
                select_2.append($("<option>").attr('value', i).text(value.name));
            }
        } else {
            $('#barplot-cluster-select-view').css('display', 'none');
            if (option == "hist") {
            } else {
                update_barplot(currently_selected_cluster);
            }
        }
    });
    $('#num-genes').change(function() {
        update_barplot(currently_selected_cluster);
    });
    update_barplot(currently_selected_cluster);

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

