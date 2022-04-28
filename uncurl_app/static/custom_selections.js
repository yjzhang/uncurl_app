//
// add a new custom colormap
function add_custom_colormap() {
    var colormap_name = window.prompt('Enter name of colormap', '');
    if (!colormap_name) {
        return 0;
    }
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
            var select = $('#label_select');
            select.empty();
            select.append('<option value="create_new_label">New label</option>');
        }
    });
}

// set the form to the given list of criteria
function set_criteria(criteria, label_name) {
    console.log('set_criteria');
    console.log(criteria);
    console.log(label_name);
    // if label_name is not an option in label_select, add it.
    var selection = $('#label_select option[value=\"'+label_name+'\"]');
    if (selection.length == 0) {
        // add label to selection
        $('#label_select').append('<option value="'+label_name+'">');
    }
    $('#label_select').val(label_name);
    $('#label_name').val(label_name);
    var form = $('#all_criteria');
    var criteria_element = $(criterion_template)[0];
    if (criteria.length == 0) {
        console.log('Criteria has length 0 in set_criteria');
        form.empty();
        var ce1 = criteria_element.cloneNode(true);
        // TODO: set label name?
        form.append(ce1);
        return 0;
    }
    form.empty();
    for (var i = 0; i < criteria.length; i++) {
        var c = criteria[i];
        console.log(c);
        var ce1 = criteria_element.cloneNode(true);
        $(ce1).find('#selection_type-1 option[value='+c.selection_type+']').attr('selected', 'selected');
        // comparisons might not be present
        if (c.selection_type == 'gene') {
            $(ce1).find('#selection_comparison-1').empty();
            $(ce1).find('#selection_comparison-1').append('<option value=">">&gt;</option>');
            $(ce1).find('#selection_comparison-1').append('<option value="<">&lt;</option>');
        }
        $(ce1).find('#selection_comparison-1 option[value="' + c.comparison +'"]').attr('selected', 'selected');
        $(ce1).find('#selection_target-1').attr('value', c.target);
        $(ce1).find('#selection_target-1').attr('list', c.selection_type);
        if (c.value) {
            $(ce1).find('#selection_value-1').toggle(true);
            $(ce1).find('#selection_value-1').attr('value', c.value);
        }
        var template = ce1.outerHTML;
        template = template.replace(/-1\"/g, '-'+String(i+1)+'"');
        template = template.replace(/\(1\)/g, '('+String(i+1)+')');
        // set selection_and_or-1
        if (i > 0) {
            template = template.replace(/\"or\"/g, '"' + c.and_or + '"');
            template = template.replace(/\"and\"/g, '"' + c.and_or + '"');
        }
        form.append(template);
    }
};

// get json colormap for the selected name 
function get_custom_colormap() {
    var colormap_name = $('#cell-color').val();
    $.ajax({url: window.location.pathname + "/get_colormap_label_criteria",
        method: 'POST',
        data: {
            name: colormap_name,
        },
    }).done(function(data) {
        if (data.startsWith('Error')) {
            $('#custom_color_map_area').toggle(false);
        } else {
            $('#custom_color_map_area').toggle(true);
            var select = $('#label_select');
            select.empty();
            data = JSON.parse(data);
            console.log(data);
            var labels = data.labels;
            // add labels to label_select
            for (var i = 0; i < labels.length; i++) {
                var l = labels[i];
                select.append('<option value="'+l.name+'">'+l.name+'</option>');
                if (i == 0) {
                    select.val(l.name);
                }
            }
            select.append('<option value="create_new_label">New label</option>');
            // set criteria to label 1
            if (data.labels.length > 0) {
                var criteria = data.labels[0].criteria;
                var label_name = data.labels[0].name;
                set_criteria(criteria, label_name);
                $('#label_color').val(data.labels[0]['color']);
            } else {
                // TODO: have some default criteria?
                set_criteria([], '');
            }
        }
    });
}

// TODO: this needs to submit a query to the server
function delete_custom_label() {
}

function delete_criterion(criterion_id) {
    $('#custom_selection_criterion-'+String(criterion_id)).remove();
}

// add a new criterion for a label
function add_custom_criterion(and_or) {
    // find number of criteria the label has...
    var num_criteria = $('.custom_selection_criterion').length;
    var id = num_criteria + 1;
    var template = criterion_template;
    template = template.replace(/-1\"/g, '-'+id+'"');
    template = template.replace(/\(1\)/g, '('+id+')');
    // set selection_and_or-1
    template = template.replace(/\"or\"/g, '"' + and_or + '"');
    template = template.replace(/\"and\"/g, '"' + and_or + '"');
    $('#all_criteria').append(template);
}

// submit label updates to server
function submit_label() {
    var colormap_name = $('#cell-color').val();
    var label_name = $('#label_name').val();
    var label_color = $('#label_color').val();
    console.log(label_color);
    console.log(label_color == '#000000');
    var criteria_form = $('#all_criteria_form').serializeArray();
    var criteria_data = {};
    $(criteria_form).each(function(index, obj){
        criteria_data[obj.name] = obj.value;
    });
    // should the array processing be done server-side or client-side?
    // whatevs just let the server handle it
    $.ajax({url: window.location.pathname + '/update_colormap_label_criteria',
        method: 'POST',
        data: {
            name: colormap_name,
            label: label_name,
            criteria: JSON.stringify(criteria_data),
            color: label_color,
        }
    }).done(function(data) {
        if (data.startsWith('Error')) {
            $("#update-area").empty();
            $("#update-area").append(data);
            return false;
        } else {
            // clear scatterplots cache, update scatterplot
            var label_data = JSON.parse(data);
            //set_criteria(label_data.criteria, label_name);
            cache.scatterplots = {};
            update_scatterplot();
        }
    });
}

// either create a new label, or set the label's criteria list to the correct list for that label.
function update_custom_label() {
    var colormap_name = $('#cell-color').val();
    var selection_val = $('#label_select').val();
    if (selection_val == 'create_new_label') {
        var label_name = window.prompt('Name for new label:');
        if (!label_name) {
            return 0;
        }
        // create a new option for label name
        $.ajax({url: window.location.pathname + "/update_colormap_label_criteria",
            method: 'POST',
            data: {
                name: colormap_name,
                label: label_name,
            },
        }).done(function(data) {
            if (data.startsWith('Error')) {
                console.log(data);
            } else {
                // add label to label_select
                $('#label_select').append('<option value="'+label_name+'">'+label_name+'</option>');
                $('#label_select').val(label_name);
                var label_data = JSON.parse(data);
                console.log(label_data);
                set_criteria(label_data.criteria, label_name);
                // add New Label option if none is available
                if ($('#label_select').children().length == 1) {
                    select.append('<option value="create_new_label">New label</option>');
                }
            }
        });
    } else {
        var label_name = selection_val;
        // get criteria for label
        $.ajax({url: window.location.pathname + "/get_colormap_label_criteria",
            method: 'POST',
            data: {
                name: colormap_name,
                label: label_name,
            },
        }).done(function(data) {
            if (data.startsWith('Error')) {
                console.log(data);
            } else {
                // add label to label_select
                var label_data = JSON.parse(data);
                set_criteria(label_data.criteria, label_name);
                console.log('label_data');
                console.log(label_data);
                $('#label_color').val(label_data['color']);
            }
        });
    }
}

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
        if ($('datalist#cluster').length == 0) {
            //var dl = $('<datalist id="cluster_datalist"></datalist>');
            // TODO: this doesn't work since the data is...
            //var n_clusters = current_scatterplot_data.length;
            //for (var i = 0; i < n_clusters; i++) {
            //    dl.append('<option value="'+i+'">');
            //}
            //dl.appendTo('body');
        }
        $('#selection_target-'+ criterion_id).attr('list', 'cluster');
        // delete gene input text box
        $('#selection_value-'+criterion_id).toggle(false);
    } else if (selection_type == 'gene') {
        // create a gene input text box after the selection
        $('#selection_value-'+criterion_id).toggle(true);
        $('#selection_value-'+criterion_id).attr('list', 'gene_names');
        comparison.append('<option value=">">&gt;</option>');
        comparison.append('<option value="<">&lt;</option>');
    } else if (selection_type == 'read_counts') {
        comparison.append('<option value=">">&gt;</option>');
        comparison.append('<option value="<">&lt;</option>');
        $('#selection_value-'+criterion_id).toggle(false);
    } else if (selection_type == 'selection') {
        comparison.append('<option value="=">=</option>');
        //$('#selection_target-'+criterion_id).attr('disabled', 'disabled');
        $('#selection_target-'+criterion_id).val(String(current_selected_cells));
        $('#selection_value-'+criterion_id).toggle(false);
    } else { // custom label
        comparison.append('<option value="=">=</option>');
        comparison.append('<option value="!=">!=</option>');
        // create a datalist for the custom label set
        if ($('datalist#'+selection_type).length == 0) {
            get_colormap_values(selection_type);
        }
        $('#selection_target-'+ criterion_id).attr('list', selection_type);
        $('#selection_value-'+criterion_id).toggle(false);
    }
}

// for all criteria, if selection_type is 'selection', then the selection_target
// is set to the current selected cells.
function set_selection_target_to_selected_cells() {
    console.log('set_selection_target_to_selected_cells');
    var num_criteria = $('.custom_selection_criterion').length;
    for (var i = 1; i <= num_criteria; i++) {
        var selection_type = $('#selection_type-'+ i).val();
        if (selection_type == 'selection') {
            $('#selection_target-'+i).val(String(current_selected_cells));
        }
    }
}

// gets available options for the colormap
function get_colormap_values(colormap) {
    var output = [];
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
            var dl = $('<datalist id="'+colormap+'"></datalist>');
            console.log('datalist length: ' + output.length);
            for (var i = 0; i < output.length; i++) {
                dl.append('<option value="'+output[i]+'">');
            }
            dl.appendTo('body');
        }
    });
}
