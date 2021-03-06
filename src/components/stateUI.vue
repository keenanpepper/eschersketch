<template>
  <div id="stateUI">
    <template v-if="params.fullUI">
      <div id="logo" class="Aligner" style="margin-top:-10px;">
        <div class="Aligner-item">
          <span class="eslogotext" style="font-variant:small-caps;margin-right:10px;">escher</span><br>
        </div>
        <div class="Aligner-item">
          <img src="/static/svg/eslogo2.svg" height="24px" style="margin-top:4px;"/>
        </div>
        <div class="Aligner-item">
          <span class="eslogotext" style="font-variant:small-caps;margin-left:10px;">sketch</span>
        </div>
        <div class="Aligner-item">
          <div class="button" :class="{selected: !params.fullUI}" style="margin-left:80%;"
              @click="toggleUI" @mouseover="setHint" hint="minimized/mobile UI mode">
            <span class="icon-shrink2"></span>
          </div>
        </div>
      </div>
    </template>
    <template v-else>
        <span class="eslogo"><img src="/static/svg/eslogo2.svg" height="20px" style="margin-bottom:-4px; padding:0px"/></span>
    </template>

    <div class="button" :class="{pulsate: pageJustLoaded}" @mousedown="help"
         key="stateui-help-button" @mouseover="setHint" hint="show help screen">
      <span class="icon-question-circle"></span>
      <!--<b>?</b>-->
    </div>

    <template v-if="params.fullUI">
      <div class="button" @mousedown="config"
           @mouseover="setHint" hint="settings e.g. turn off these hints">
        <span class="icon-cog"></span>
      </div>
    </template>

    <template v-if="!params.fullUI">
      <div class="button" :class="{selected: params.showTool}" @mousedown="toggleTool"
            key="stateui-tool-button"
           @mouseover="setHint" hint="show drawing tools">
        <span class="icon-quill"></span>
      </div>
      <div class="button" :class="{selected: params.showColor}" @mousedown="toggleColor"
           key="stateui-color-button"
           @mouseover="setHint" hint="show color selection">
        <span class="icon-palette"></span>
      </div>
      <div class="button" :class="{selected: params.showSymm}" @mousedown="toggleSymm"
           key="stateui-symm-button"
           @mouseover="setHint" hint="show symmetry selection">
        <span class="icon-symmetries"></span>
      </div>
      <div class="button" :class="{selected: params.showFile}" @mousedown="toggleFile"
           key="stateui-file-button"
           @mouseover="setHint" hint="show save and share options">
        <span class="icon-folder-download"></span>
      </div>
    </template>

    <div class="button" @mousedown="undo" key="stateui-undo"
         @mouseover="setHint" hint="step back one drawing operation">
      <span class="icon-undo"></span>
    </div>
    <div class="button" @mousedown="redo" key="stateui-redo"
         @mouseover="setHint" hint="re-step forward one drawing operation">
      <span class="icon-redo"></span>
    </div>

    <div class="button" :class="{armed: armed}" @mousedown="reset" key="stateui-reset"
         @mouseover="setHint" hint="reset drawing (double click to confirm)">
      <template v-if="armed"><span class="icon-bin"></span>?</template>
      <template v-else><span class="icon-bin"></span></template>
    </div>


    <template v-if="!params.fullUI">
      <div class="button" @click="toggleUI" key="stateui-enlarge-button"
           @mouseover="setHint" hint="full UI mode">
        <span class="icon-enlarge2"></span>
      </div>
    </template>

</div>
</template>

<script>
import {gS} from '../main.js';

export default {
  props: ['params'],
  data: function(){ return {toggled: false, armed: false, pageJustLoaded:true}; },
  components: {},
  mounted: function() {
    setTimeout(() => this.pageJustLoaded=false, 10000);
  },
  computed:{
    toggleClass: function() {
      if(this.toggled) {
        return "alarm"
      }
    },
    whichButton: function(name){
      if(name=="tool"){
        return {}
      }
    }
  },
  methods: {
    undo: function(){ gS.$emit('undo'); },
    redo: function(){ gS.$emit('redo'); },
    reset: function(){
      if(this.armed){
        gS.$emit('reset');
        this.armed=false;
      } else {
        this.armed=true;
        setTimeout(() => this.armed=false, 1000);
      }
    },
    toggleUI: function(){ gS.$emit('toggleUI'); },
    toggleTool: function(){
      gS.$emit('toggleParam', 'showTool');
      gS.$emit('toggleParam', 'showLine');
    },
    toggleColor: function(){
      gS.$emit('toggleParam', 'showColor');
    },
    toggleSymm: function(){
      gS.$emit('toggleParam', 'showSymm');
    },
    toggleFile: function(){
      gS.$emit('toggleParam', 'showFile');
      gS.$emit('toggleParam', 'showNetwork');
    },
    help: function(){ gS.$emit('help'); },
    config: function(){ gS.$emit('config'); },
    setHint: function(e) {
      if(e.target.attributes.hint){
        gS.$emit('setHint', e.target.attributes.hint.value);
      }
    },
  },
}
</script>

<style scoped>

[tooltip]:before {
    position : absolute;
    background-color: #fff;
    content : attr(tooltip);
    opacity : 0;
}
[tooltip]:hover:before {
    opacity : 1;
}

.flex-container {
    height: 100%;
    padding: 0;
    margin: 0;
    display: -webkit-box;
    display: -moz-box;
    display: -ms-flexbox;
    display: -webkit-flex;
    display: flex;
    align-items: center;
    justify-content: center;
}
.row {
    width: auto;
    border: 1px solid blue;
}
.flex-item {
    background-color: tomato;
    padding: 5px;
    width: 20px;
    height: 20px;
    margin: 10px;
    line-height: 20px;
    color: white;
    font-weight: bold;
    font-size: 2em;
    text-align: center;
}

/* flexbox v-center */
.Aligner {
  display: flex;
  align-items: center;
  justify-content: center;
}
.Aligner-item {
  max-width: 50%;
}
.Aligner-item--top {
  align-self: flex-start;
}
.Aligner-item--bottom {
  align-self: flex-end;
}



@-webkit-keyframes color_change {
	from { background-color: #eeeeee; }
	to { background-color: #ff8888; }
}
@-moz-keyframes color_change {
	from { background-color: #eeeeee; }
	to { background-color: #ff8888; }
}
@-ms-keyframes color_change {
	from { background-color: #eeeeee; }
	to { background-color: #ff8888; }
}
@-o-keyframes color_change {
	from { background-color: #eeeeee; }
	to { background-color: #ff8888; }
}
@keyframes color_change {
	from { background-color: #eeeeee; }
	to { background-color: #ff8888; }
}
.pulsate {
	-webkit-animation: color_change 0.5s 20 alternate;
	-moz-animation: color_change 0.5s 20 alternate;
	-ms-animation: color_change 0.5s 20 alternate;
	-o-animation: color_change 0.5s 20 alternate;
	animation: color_change 0.5s 20 alternate;
}

</style>
