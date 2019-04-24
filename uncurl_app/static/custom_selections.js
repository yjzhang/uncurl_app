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
        }
    });
};

// set the form to the given list of criteria
function set_criteria(criteria) {
    var form = $('#all_criteria');
    var criteria_element = $(criterion_template)[0];
    if (criteria.length == 0) {
        console.log('Error: criteria has length 0 in set_criteria');
        return 0;
    }
    form.empty();
    for (var i = 0; i < criteria.length; i++) {
        var c = criteria[i];
        console.log(c);
        var ce1 = criteria_element.cloneNode(true);
        $(ce1).find('#selection_type-1 option[value='+c.selection_type+']').attr('selected', 'selected');
        $(ce1).find('#selection_comparison-1 option[value=\\' + c.comparison +']').attr('selected', 'selected');
        $(ce1).find('#selection_target-1').attr('value', c.target);
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
            select.append('<option value="new_label">New label</option>');
            // set criteria to label 1
            if (data.labels.length > 0) {
                var criteria = data.labels[0].criteria;
                set_criteria(criteria);
            } else {
                // TODO: have some default criteria?
            }
        }
    });
};

// TODO
function delete_custom_label() {
};

// TODO
function delete_criterion(criterion_id) {
};

// add a new criterion for a label
function add_custom_criterion(and_or) {
    // find number of criteria the label has...
    var num_criteria = $('.custom_selection_criterion').length;
    var id = num_criteria + 1;
    var template = $('#custom_selection_criterion-1')[0].outerHTML;
    template = template.replace(/-1\"/g, '-'+id+'"');
    template = template.replace(/\(1\)/g, '('+id+')');
    // set selection_and_or-1
    template = template.replace(/\"or\"/g, '"' + and_or + '"');
    template = template.replace(/\"and\"/g, '"' + and_or + '"');
    console.log(template);
    $('#all_criteria').append(template);
};

// TODO: submit label updates to server
function submit_label() {
    var colormap_name = $('#cell-color').val();
    var selection_val = $('#label_select').val();
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
            label: selection_val,
            criteria: JSON.stringify(criteria_data),
        }
    }).done(function(data) {
        if (data.startsWith('Error')) {
            $("#update-area").empty();
            $("#update-area").append(data);
            return false;
        } else {
            // TODO: clear scatterplots cache, update scatterplot
            cache.scatterplots = {};
            update_scatterplot();
        }
    });
};

// either create a new label, or set the label's criteria list to the correct list for that label.
function update_custom_label() {
    var colormap_name = $('#cell-color').val();
    var selection_val = $('#label_select').val();
    if (selection_val == 'new_label') {
        var label_name = window.prompt('Name for new label:');
        if (!label_name) {
            return 0;
        }
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
                var label_data = JSON.parse(data);
                set_criteria(label_data);
            }
        });
    } else {
        // get criteria for label
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
                var label_data = JSON.parse(data);
                set_criteria(label_data);
            }
        });
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
