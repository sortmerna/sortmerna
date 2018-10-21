/**
 * FILE: main.js
 * Created: Oct 08, 2018 Mon
 */

LIGHT_STYLE = "./css/light.css"
LIGHT       = "Light"
DARK_STYLE  = "./css/dark.css"
DARK        = "Dark"

function setTheme() {
    console.log(this.options[this.selectedIndex].text);
    if (this.selectedIndex == 1) {
        document.getElementById("color_css").href = DARK_STYLE
    }
    else if (this.selectedIndex == 0) {
        document.getElementById("color_css").href = LIGHT_STYLE
    }
 }

function buildMenu() {
    var sel = document.createElement("select");
    sel.setAttribute("id", "menu_select");
    sel.addEventListener("change", setTheme)
    document.getElementById("menu").appendChild(sel);
    // Light
    opt = document.createElement("option");
    opt.text = LIGHT
    sel.add(opt)
    // Dark
    var opt = document.createElement("option");
    opt.text = DARK
    sel.add(opt)
}

buildMenu();