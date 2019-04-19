

// add a new custom colormap
function add_custom_colormap() {
    var colormap_name = window.prompt('Enter name of colormap', '');
    $.ajax({url: window.location.pathname + "/custom_color_map",
        method: 'POST',
        data: {
            name: colormap_name,
        },
    }).done(function(data) {
        if (data.startsWith('Error')) {
            $('#update-area').empty();
            $('#update-area').append(data);
        } else {
            // add label to cell-color
            $('#cell-color').append('<option value="' + colormap_name + '">' + colormap_name + '</option>');
            // set name to selected name
            $('#cell-color').val(colormap_name);
            $('#custom_color_map_area').toggle(true);
        }
    });
};

// TODO: set form criteria for the specified criteria
// criteria is a list of criteria...
function set_criteria(criteria) {
    for (var i = 0; i < criteria.length; i++) {
        var c = criteria[i];
    }
};

// TODO: get json colormap for the selected name 
function get_custom_colormap() {
    var colormap_name = $('#cell-color').val();
    $.ajax({url: window.location.pathname + "/get_colormap_label_criteria",
        method: 'POST',
        data: {
            name: colormap_name,
        },
    }).done(function(data) {
        if (data.startsWith('Error')) {
        } else {
            $('#label_select').empty();
            data = JSON.parse(data);
            var labels = data.labels;
            // add labels to label_select
            for (var i = 0; i < labels.length; i++) {
                var l = labels[i];
                $('#label_select').append('<option value="'+l.name+'">'+l.name+'</option>');
                if (i == 0) {
                    $('#label_select').val(l.name);
                }
            }
            $('#label_select').append('<option value="new_label">New label</option>');
            // TODO: set criteria to label 1
            var criteria = data.labels[0].criteria;
            set_criteria(criteria);
        }
    });
};

// TODO
function delete_custom_label() {
};


// add a new criterion for a label
function add_custom_criterion(and_or) {
    // find number of criteria the label has...
    var num_criteria = $('.custom_selection_criterion').length;
    var id = num_criteria + 1;
    var template = $('#custom_selection_criterion-1').html();
    template.replace(/-1\"/g, '-'+id+'"');
    template.replace(/update_custom_criterion\(1\)/g, 'update_custom_criterion('+id+')');
    // set selection_and_or-1
    template.replace(/\"or\"/g, '"' + and_or + '"');
    $('#all_criteria').append(template);
};

// TODO: submit label updates to server
function submit_label() {
    var colormap_name = $('#cell-color').val();
    var selection_val = $('#label_select').val();
    var criteria_form = $('#all_criteria_form').serializeArray();
};

// TODO: delete a criterion for a label
function delete_custom_criterion(label_id, criterion_id) {
};

// TODO: either create a new label, or set the label's criteria list to the correct list for that label.
function update_custom_label() {
    var colormap_name = $('#cell-color').val();
    var selection_val = $('#label_select').val();
    if (selection_val == 'new_label') {
        var label_name = window.prompt('Name for new label:');
        $.ajax({url: window.location.pathname + "/update_colormap_label_criteria",
            method: 'POST',
            data: {
                name: colormap_name,
                label: label_name,
            },
        }).done(function(data) {
            if (data.startsWith('Error')) {
            } else {
                // add label to label_select
                $('#label_select').append('<option value="'+label_name+'">'+label_name+'</option>');
                $('#label_select').val(label_name);
                // TODO: set criteria
                var label_data = JSON.parse(data);
            }
        });
    } else {
        // TODO: get criteria for label
    }
};

// change the available options
function update_custom_criterion(criterion_id) {
    var selection_type = $('#selection_type-'+ criterion_id).val();
    var comparison = $('#selection_comparison-'+ criterion_id);
    comparison.empty();
    // if cluster, the comparison can be 
    if (selection_type == 'cluster') {
        comparison.append('<option value="=">=</option>');
        comparison.append('<option value="!=">!=</option>');
        // create a datalist for all clusters
        if ($('datalist#cluster_datalist').length == 0) {
            var dl = $('<datalist id="cluster_datalist"></datalist>');
            var n_clusters = current_scatterplot_data.length;
            for (var i = 0; i < n_clusters; i++) {
                dl.append('<option value="'+i+'">');
            }
            dl.appendTo('body');
        }
        $('#selection_target-'+ criterion_id).attr('datalist', 'cluster_datalist');
    } else if (selection_type == 'gene') {
        // TODO: create a gene input
        comparison.append('<option value=">=">greater than</option>');
        comparison.append('<option value="!=">less than</option>');
    } else if (selection_type == 'read_counts') {
        comparison.append('<option value=">=">greater than</option>');
        comparison.append('<option value="!=">less than</option>');
    } else { // custom label
        comparison.append('<option value="=">=</option>');
        comparison.append('<option value="!=">!=</option>');
        // create a datalist for the custom label set
        if ($('datalist#'+selection_type).length == 0) {
            var options = get_colormap_values(selection_type);
            var dl = $('<datalist id="'+selection_type+'"></datalist>');
            for (var i = 0; i < options.length; i++) {
                dl.append('<option value="'+options[i]+'">');
            }
            dl.appendTo('body');
        }
        $('#selection_target-'+ criterion_id).attr('datalist', selection_type);
    }
};

// gets available options for the colormap
function get_colormap_values(colormap) {
    var output = {};
    $.ajax({url: window.location.pathname + "/get_colormap_values",
        method: 'POST',
        data: {
            name: colormap,
        },
    }).done(function(data) {
        if (data.startsWith('Error')) {
            $('#update-area').empty();
            $('#update-area').append(data);
        } else {
            output = JSON.parse(data);
            console.log(output);
        }
    });
    console.log(output);
    return output;
};
