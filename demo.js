'use strict';
console.log('start demo');

var canvas = document.querySelector('#canvas');
var polyline = canvas.querySelector('#polyline');
var path = canvas.querySelector('#path');
var isPathclose = document.querySelector('#is-pathclose');
var points = [];
var isDrag = false;

function drawPolyline(points) {
  var points_ = [];
  points.forEach(function(point) {
    points_.push(point.join(',')); 
  });
  polyline.setAttribute('points', points_.join(' '));
}

function drawPath(segments) {
  var d = [];
  segments.forEach(function(segment, index) {
    if (index === 0) {
      d.push('M' + segment[0].x + ',' + segment[0].y + ' C');
    }
    segment.forEach(function(point, index_) {
      if (index_ === 0) {
        return;
      }
      d.push(point.x + ',' + point.y);
    });
  });
  if (segments.length > 0 && isPathclose.checked) {
    d.push('z');
  }
  d = d.join(' ');
  console.log('d="' + d + '"');
  path.setAttribute('d', d);
}

function addPoints(x, y) {
  points.push([x, y]);
  drawPolyline(points);
}

function onMouseDown() {
  points = [];
  isDrag = true;
  drawPolyline(points);
  drawPath([]);
}

function onMouseMove(e) {
  if (isDrag) {
    var pointer = typeof e.changedTouches !== 'undefined' ?
      e.changedTouches[0] : e;
    var x = pointer.clientX; 
    var y = pointer.clientY; 
    mouseLogging(x, y);
    addPoints(x, y);
    e.preventDefault();
  }
}

function onMouseUp() {
  if (isDrag) {
    isDrag = false;
    console.log('points: ' + JSON.stringify(points));
    if (isPathclose.checked) {
      points.push(points[0]);
    }
    var segments = polyline2bezier(points, 10);
    drawPath(segments);
    points = [];
    drawPolyline(points);
  }
}

// add event
canvas.addEventListener('mousedown', onMouseDown, false);
document.addEventListener('mousemove', onMouseMove, false);
document.addEventListener('mouseup', onMouseUp, false);
canvas.addEventListener('touchstart', onMouseDown, false);
document.addEventListener('touchmove', onMouseMove, false);
document.addEventListener('touchend', onMouseUp, false);

function mouseLogging(x, y) {
  console.log({x: x, y: y});
}
