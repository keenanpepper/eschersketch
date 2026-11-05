//------------------------------------------------------------------------------
//
// Eschersketch - A drawing program for exploring symmetrical designs
//
//
// Copyright (c) 2026 Anselm Levskaya (http://anselmlevskaya.com)
// Licensed under the MIT (http://www.opensource.org/licenses/mit-license.php)
// license.
//
//------------------------------------------------------------------------------

// DRAWING GLOBALS
import {gS,
        livecanvas, lctx, canvas, ctx,
        drawKeyToOrderMap, symmetrySet,
        startBrush, loadBrush, clearBrush, liveStyle, liveSymmetry,
        commitOp
       } from './main';
import _ from 'underscore';
import {l2dist} from './math_utils';

import {drawHitCircle} from './canvas_utils';
import {mappedPath} from './mappedPath';

// Draw Single Line Segments
//------------------------------------------------------------------------------
export class LineOp {
    constructor(symmState, ctxStyle, points) {
    this.tool = "line";
    this.points = points;
    this.ctxStyle = ctxStyle;
    this.symmState = symmState;
  }

  render(ctx){
    _.assign(ctx, this.ctxStyle);
    //gS.$emit('symmUpdate', this.symmState);
    for (let af of symmetrySet(this.symmState)) {
      const path = mappedPath(ctx, af);
      ctx.beginPath();
      path.moveTo(this.points[0]);
      path.lineTo(this.points[1]);
      ctx.stroke();
    }
  }
}

//State Labels
const _INIT_ = 0;
const _OFF_  = 1;
const _ON_   = 2;
const _MOVESTART_ = 3;
const _MOVEEND_ = 4;

export class LineTool {
  constructor() {
    this.points = [[0,0],[0,0]];
    this.state = _INIT_;
    this.hitRadius = 4;
    this.actions = [
      {name: "cancel", desc: "cancel line",    icon: "icon-cross",     key: "Escape"},
      {name: "commit", desc: "start new (automatic on new click)", icon: "icon-checkmark", key: "Enter"},
    ];
  }

  liverender() {
    _.assign(lctx, liveStyle(this));
    lctx.clearRect(0, 0, canvas.width, canvas.height);
    for (let af of liveSymmetry(this)) {
      const path = mappedPath(lctx, af);
      lctx.beginPath();
      path.moveTo(this.points[0]);
      path.lineTo(this.points[1]);
      lctx.stroke();
    }
    drawHitCircle(lctx, this.points[0][0]-0.5, this.points[0][1]-0.5, this.hitRadius);
    drawHitCircle(lctx, this.points[1][0]-0.5, this.points[1][1]-0.5, this.hitRadius);
  }

  enter(op){
    if(op){
        loadBrush(this, op);
        this.points = op.points;
        this.state = _OFF_;
        this.liverender();
    } else{
      this.points = [[0,0],[0,0]];
      this.state = _INIT_;
    }
  }

  exit(){
    clearBrush(this);
      this.points = [[0,0],[0,0]];
      this.state = _INIT_;
  }

  commit() {
    if(this.state == _INIT_){return;}
    let ctxStyle = _.clone(liveStyle(this));
    let symmState = _.clone(this.symmState || gS.symmState);
    commitOp(new LineOp(symmState, ctxStyle, this.points));
    clearBrush(this);
    lctx.clearRect(0, 0, livecanvas.width, livecanvas.height);
    this.points = [[0,0],[0,0]];
    this.state = _INIT_;
  }

  cancel() {
    clearBrush(this);
    lctx.clearRect(0, 0, livecanvas.width, livecanvas.height);
    this.state = _INIT_;
    this.points = [[0,0],[0,0]];
  }

  mouseDown(e) {
    let rect = livecanvas.getBoundingClientRect();
    let pt = [e.clientX-rect.left, e.clientY-rect.top];
    if(l2dist(pt, this.points[0])<this.hitRadius) {
      this.state = _MOVESTART_;
    } else if(l2dist(pt, this.points[1])<this.hitRadius) {
      this.state = _MOVEEND_;
    } else {
      if(this.state==_OFF_) {
        this.commit();
      }
      this.state = _ON_;
      startBrush(this);
      this.points = [pt, pt];
    }
  }

  mouseMove(e) {
    let rect = livecanvas.getBoundingClientRect();
    let pt = [e.clientX-rect.left, e.clientY-rect.top];
    if (this.state == _ON_) {
        this.points[1] = pt;
        this.liverender();
    }
    else if (this.state == _MOVESTART_) {
      this.points[0] = pt;
      this.liverender();
    }
    else if (this.state == _MOVEEND_) {
      this.points[1] = pt;
      this.liverender();
    }
  }

  mouseUp(e) {
    if(this.state===_INIT_){return;} //edgecase of accidental mouseup before drawing
    this.state = _OFF_;
  }

  keyDown(e) {
    if(e.target.type){return;} // don't interfere with input UI key-events
    for(let action of this.actions){
      if(action.key == e.code){
        this[action.name]();
      }
    }
  }

}
