alert('If you see this alert, then your custom JavaScript script has run!')
$(document).ready(function() {
    alert("you will only see this after the file is loaded and jquery is imported")
    $( ".sidebar" ).resizable({
        handles: "e, w"
    });
    $( ".graph" ).resizable({
        handles: "e, w"
    });
    $( "summary" ).click(function () {
        $( this ).toggleClass("toggle")
    })
}
)
