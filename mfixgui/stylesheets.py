DEFAULT_STYLESHEET = """

QSplitter::handle:horizontal {
    image: url(:images/left_right_arrow.svg);
    background-color: #e0e0e0;
    width: 6px;
}
QSplitter::handle:vertical {
    image: url(:images/up_down_arrow.svg);
    background-color: #e0e0e0;
    height: 6px;
}
"""



OSX_STYLESHEET = DEFAULT_STYLESHEET + """

QToolButton:checked {
  border: 2px solid #8f8f91;
  background-color: rgb(175,175,175);
  border-radius: 4px
}
QToolButton {
  background-color: transparent;
}
QToolButton:hover {
  background-color: rgb(200,200,200);
  border-radius: 4px
}
"""
