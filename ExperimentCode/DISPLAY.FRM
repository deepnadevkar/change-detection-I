VERSION 5.00
Object = "{C30E0AD3-B1BA-11CE-ABC6-F5B2E79D9E3F}#6.0#0"; "lead32.ocx"
Begin VB.Form frmDisplay 
   Appearance      =   0  'Flat
   BackColor       =   &H00000000&
   BorderStyle     =   0  'None
   ClientHeight    =   4410
   ClientLeft      =   780
   ClientTop       =   2220
   ClientWidth     =   7005
   ControlBox      =   0   'False
   KeyPreview      =   -1  'True
   LinkTopic       =   "Form1"
   MaxButton       =   0   'False
   MinButton       =   0   'False
   PaletteMode     =   1  'UseZOrder
   ScaleHeight     =   4410
   ScaleMode       =   0  'User
   ScaleWidth      =   7334.832
   ShowInTaskbar   =   0   'False
   Visible         =   0   'False
   Begin VB.TextBox Text1 
      Height          =   495
      Left            =   2880
      TabIndex        =   4
      Text            =   "Text1"
      Top             =   3840
      Visible         =   0   'False
      Width           =   1215
   End
   Begin LeadLib.Lead leadStimuli 
      Height          =   492
      Left            =   0
      TabIndex        =   0
      Top             =   0
      Visible         =   0   'False
      Width           =   480
      _Version        =   393216
      _ExtentX        =   843
      _ExtentY        =   882
      _StockProps     =   97
      BackColor       =   16711680
      AutoScroll      =   0   'False
      ScaleHeight     =   495
      ScaleMode       =   1
      ScaleWidth      =   480
      DataField       =   ""
   End
   Begin LeadLib.Lead leadProbe 
      Height          =   495
      Left            =   1440
      TabIndex        =   1
      Top             =   3600
      Visible         =   0   'False
      Width           =   480
      _Version        =   393216
      _ExtentX        =   843
      _ExtentY        =   882
      _StockProps     =   97
      BackColor       =   255
      AutoScroll      =   0   'False
      ScaleHeight     =   495
      ScaleMode       =   1
      ScaleWidth      =   480
      DataField       =   ""
   End
   Begin LeadLib.Lead Lead1 
      Height          =   495
      Left            =   480
      TabIndex        =   6
      Top             =   0
      Visible         =   0   'False
      Width           =   480
      _Version        =   393216
      _ExtentX        =   843
      _ExtentY        =   882
      _StockProps     =   97
      BackColor       =   65535
      AutoScroll      =   0   'False
      ScaleHeight     =   495
      ScaleMode       =   1
      ScaleWidth      =   480
      DataField       =   ""
   End
   Begin LeadLib.Lead Lead2 
      Height          =   495
      Left            =   960
      TabIndex        =   7
      Top             =   0
      Visible         =   0   'False
      Width           =   480
      _Version        =   393216
      _ExtentX        =   843
      _ExtentY        =   882
      _StockProps     =   97
      BackColor       =   16744576
      AutoScroll      =   0   'False
      ScaleHeight     =   495
      ScaleMode       =   1
      ScaleWidth      =   480
      DataField       =   ""
   End
   Begin LeadLib.Lead Lead3 
      Height          =   495
      Left            =   1440
      TabIndex        =   8
      Top             =   0
      Visible         =   0   'False
      Width           =   480
      _Version        =   393216
      _ExtentX        =   843
      _ExtentY        =   882
      _StockProps     =   97
      BackColor       =   16776960
      AutoScroll      =   0   'False
      ScaleHeight     =   495
      ScaleMode       =   1
      ScaleWidth      =   480
      DataField       =   ""
   End
   Begin LeadLib.Lead Lead4 
      Height          =   495
      Left            =   0
      TabIndex        =   9
      Top             =   480
      Visible         =   0   'False
      Width           =   480
      _Version        =   393216
      _ExtentX        =   843
      _ExtentY        =   882
      _StockProps     =   97
      BackColor       =   33023
      AutoScroll      =   0   'False
      ScaleHeight     =   495
      ScaleMode       =   1
      ScaleWidth      =   480
      DataField       =   ""
   End
   Begin LeadLib.Lead Lead5 
      Height          =   495
      Left            =   480
      TabIndex        =   10
      Top             =   480
      Visible         =   0   'False
      Width           =   480
      _Version        =   393216
      _ExtentX        =   843
      _ExtentY        =   882
      _StockProps     =   97
      BackColor       =   12648384
      AutoScroll      =   0   'False
      ScaleHeight     =   495
      ScaleMode       =   1
      ScaleWidth      =   480
      DataField       =   ""
   End
   Begin LeadLib.Lead Lead6 
      Height          =   495
      Left            =   960
      TabIndex        =   11
      Top             =   480
      Visible         =   0   'False
      Width           =   480
      _Version        =   393216
      _ExtentX        =   843
      _ExtentY        =   882
      _StockProps     =   97
      BackColor       =   16761087
      AutoScroll      =   0   'False
      ScaleHeight     =   495
      ScaleMode       =   1
      ScaleWidth      =   480
      DataField       =   ""
   End
   Begin LeadLib.Lead Lead7 
      Height          =   492
      Left            =   1440
      TabIndex        =   12
      Top             =   480
      Visible         =   0   'False
      Width           =   480
      _Version        =   393216
      _ExtentX        =   843
      _ExtentY        =   882
      _StockProps     =   97
      BackColor       =   65280
      AutoScroll      =   0   'False
      ScaleHeight     =   495
      ScaleMode       =   1
      ScaleWidth      =   480
      DataField       =   ""
   End
   Begin LeadLib.Lead Lead8 
      Height          =   492
      Left            =   0
      TabIndex        =   13
      Top             =   960
      Visible         =   0   'False
      Width           =   480
      _Version        =   393216
      _ExtentX        =   843
      _ExtentY        =   882
      _StockProps     =   97
      BackColor       =   8421631
      AutoScroll      =   0   'False
      ScaleHeight     =   495
      ScaleMode       =   1
      ScaleWidth      =   480
      DataField       =   ""
   End
   Begin LeadLib.Lead Lead9 
      Height          =   492
      Left            =   480
      TabIndex        =   14
      Top             =   960
      Visible         =   0   'False
      Width           =   480
      _Version        =   393216
      _ExtentX        =   843
      _ExtentY        =   882
      _StockProps     =   97
      BackColor       =   8421504
      AutoScroll      =   0   'False
      ScaleHeight     =   495
      ScaleMode       =   1
      ScaleWidth      =   480
      DataField       =   ""
   End
   Begin LeadLib.Lead Lead10 
      Height          =   492
      Left            =   960
      TabIndex        =   15
      Top             =   960
      Visible         =   0   'False
      Width           =   480
      _Version        =   393216
      _ExtentX        =   843
      _ExtentY        =   882
      _StockProps     =   97
      BackColor       =   12632064
      AutoScroll      =   0   'False
      ScaleHeight     =   495
      ScaleMode       =   1
      ScaleWidth      =   480
      DataField       =   ""
   End
   Begin LeadLib.Lead Lead11 
      Height          =   492
      Left            =   1440
      TabIndex        =   16
      Top             =   960
      Visible         =   0   'False
      Width           =   480
      _Version        =   393216
      _ExtentX        =   843
      _ExtentY        =   882
      _StockProps     =   97
      BackColor       =   12583104
      AutoScroll      =   0   'False
      ScaleHeight     =   495
      ScaleMode       =   1
      ScaleWidth      =   480
      DataField       =   ""
   End
   Begin LeadLib.Lead Lead12 
      Height          =   492
      Left            =   0
      TabIndex        =   17
      Top             =   1440
      Visible         =   0   'False
      Width           =   480
      _Version        =   393216
      _ExtentX        =   843
      _ExtentY        =   882
      _StockProps     =   97
      BackColor       =   8438015
      AutoScroll      =   0   'False
      ScaleHeight     =   495
      ScaleMode       =   1
      ScaleWidth      =   480
      DataField       =   ""
   End
   Begin LeadLib.Lead Lead13 
      Height          =   492
      Left            =   480
      TabIndex        =   18
      Top             =   1440
      Visible         =   0   'False
      Width           =   480
      _Version        =   393216
      _ExtentX        =   843
      _ExtentY        =   882
      _StockProps     =   97
      BackColor       =   128
      AutoScroll      =   0   'False
      ScaleHeight     =   495
      ScaleMode       =   1
      ScaleWidth      =   480
      DataField       =   ""
   End
   Begin LeadLib.Lead Lead14 
      Height          =   492
      Left            =   960
      TabIndex        =   19
      Top             =   1440
      Visible         =   0   'False
      Width           =   480
      _Version        =   393216
      _ExtentX        =   843
      _ExtentY        =   882
      _StockProps     =   97
      BackColor       =   8421376
      AutoScroll      =   0   'False
      ScaleHeight     =   495
      ScaleMode       =   1
      ScaleWidth      =   480
      DataField       =   ""
   End
   Begin LeadLib.Lead Lead15 
      Height          =   492
      Left            =   1440
      TabIndex        =   20
      Top             =   1440
      Visible         =   0   'False
      Width           =   480
      _Version        =   393216
      _ExtentX        =   843
      _ExtentY        =   882
      _StockProps     =   97
      BackColor       =   8388608
      AutoScroll      =   0   'False
      ScaleHeight     =   495
      ScaleMode       =   1
      ScaleWidth      =   480
      DataField       =   ""
   End
   Begin VB.Label lblDiff 
      AutoSize        =   -1  'True
      BackColor       =   &H00FFFFFF&
      Height          =   435
      Left            =   3360
      TabIndex        =   2
      Top             =   2520
      Visible         =   0   'False
      Width           =   645
   End
   Begin VB.Label lblDiffTouch 
      BackColor       =   &H00000000&
      Height          =   1395
      Left            =   4560
      TabIndex        =   5
      Top             =   840
      Visible         =   0   'False
      Width           =   1455
   End
   Begin VB.Label lblProbe 
      BackColor       =   &H00FFFFFF&
      Height          =   2715
      Left            =   4320
      TabIndex        =   3
      Top             =   240
      Visible         =   0   'False
      Width           =   2685
   End
End
Attribute VB_Name = "frmDisplay"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = False
Option Explicit
Dim l As Integer


Private Sub Form_Activate()

    Refresh
    
End Sub

Private Sub build_array()
    Set stimulus_array(1) = leadStimuli
    Set stimulus_array(2) = Lead1
    Set stimulus_array(3) = Lead2
    Set stimulus_array(4) = Lead3
    Set stimulus_array(5) = Lead4
    Set stimulus_array(6) = Lead5
    Set stimulus_array(7) = Lead6
    Set stimulus_array(8) = Lead7
    Set stimulus_array(9) = Lead8
    Set stimulus_array(10) = Lead9
    Set stimulus_array(11) = Lead10
    Set stimulus_array(12) = Lead11
    Set stimulus_array(13) = Lead12
    Set stimulus_array(14) = Lead13
    Set stimulus_array(15) = Lead14
    Set stimulus_array(16) = Lead15
End Sub

Private Sub Form_Load()
    build_array
    leadStimuli.BackColor = COLOR_BLACK
    Lead1.BackColor = COLOR_BLACK
    Lead2.BackColor = COLOR_BLACK
    Lead3.BackColor = COLOR_BLACK
    Lead4.BackColor = COLOR_BLACK
    Lead5.BackColor = COLOR_BLACK
    Lead6.BackColor = COLOR_BLACK
    Lead7.BackColor = COLOR_BLACK
    Lead8.BackColor = COLOR_BLACK
    Lead9.BackColor = COLOR_BLACK
    Lead10.BackColor = COLOR_BLACK
    Lead11.BackColor = COLOR_BLACK
    Lead12.BackColor = COLOR_BLACK
    Lead13.BackColor = COLOR_BLACK
    Lead14.BackColor = COLOR_BLACK
    Lead15.BackColor = COLOR_BLACK
    leadProbe.BackColor = COLOR_RED
    
    leadStimuli.Top = STIMULI_TOP * Screen.TwipsPerPixelY
    leadStimuli.Left = STIMULI_LEFT * Screen.TwipsPerPixelX
    leadStimuli.Height = STIMULI_HEIGHT * Screen.TwipsPerPixelY
    leadStimuli.Width = STIMULI_WIDTH * Screen.TwipsPerPixelX
    
    Lead1.Top = LEAD1_TOP * Screen.TwipsPerPixelY
    Lead1.Left = LEAD1_LEFT * Screen.TwipsPerPixelX
    Lead1.Height = LEAD1_HEIGHT * Screen.TwipsPerPixelY
    Lead1.Width = LEAD1_WIDTH * Screen.TwipsPerPixelX
    
    Lead2.Top = LEAD2_TOP * Screen.TwipsPerPixelY
    Lead2.Left = LEAD2_LEFT * Screen.TwipsPerPixelX
    Lead2.Height = LEAD2_HEIGHT * Screen.TwipsPerPixelY
    Lead2.Width = LEAD2_WIDTH * Screen.TwipsPerPixelX
    
    Lead3.Top = LEAD3_TOP * Screen.TwipsPerPixelY
    Lead3.Left = LEAD3_LEFT * Screen.TwipsPerPixelX
    Lead3.Height = LEAD3_HEIGHT * Screen.TwipsPerPixelY
    Lead3.Width = LEAD3_WIDTH * Screen.TwipsPerPixelX
    
    Lead4.Top = LEAD4_TOP * Screen.TwipsPerPixelY
    Lead4.Left = LEAD4_LEFT * Screen.TwipsPerPixelX
    Lead4.Height = LEAD4_HEIGHT * Screen.TwipsPerPixelY
    Lead4.Width = LEAD4_WIDTH * Screen.TwipsPerPixelX
    
    Lead5.Top = LEAD5_TOP * Screen.TwipsPerPixelY
    Lead5.Left = LEAD5_LEFT * Screen.TwipsPerPixelX
    Lead5.Height = LEAD5_HEIGHT * Screen.TwipsPerPixelY
    Lead5.Width = LEAD5_WIDTH * Screen.TwipsPerPixelX
    
    Lead6.Top = LEAD6_TOP * Screen.TwipsPerPixelY
    Lead6.Left = LEAD6_LEFT * Screen.TwipsPerPixelX
    Lead6.Height = LEAD6_HEIGHT * Screen.TwipsPerPixelY
    Lead6.Width = LEAD6_WIDTH * Screen.TwipsPerPixelX
    
    Lead7.Top = LEAD7_TOP * Screen.TwipsPerPixelY
    Lead7.Left = LEAD7_LEFT * Screen.TwipsPerPixelX
    Lead7.Height = LEAD7_HEIGHT * Screen.TwipsPerPixelY
    Lead7.Width = LEAD7_WIDTH * Screen.TwipsPerPixelX
    
    Lead8.Top = LEAD8_TOP * Screen.TwipsPerPixelY
    Lead8.Left = LEAD8_LEFT * Screen.TwipsPerPixelX
    Lead8.Height = LEAD8_HEIGHT * Screen.TwipsPerPixelY
    Lead8.Width = LEAD8_WIDTH * Screen.TwipsPerPixelX
    
    Lead9.Top = LEAD9_TOP * Screen.TwipsPerPixelY
    Lead9.Left = LEAD9_LEFT * Screen.TwipsPerPixelX
    Lead9.Height = LEAD9_HEIGHT * Screen.TwipsPerPixelY
    Lead9.Width = LEAD9_WIDTH * Screen.TwipsPerPixelX
    
    Lead10.Top = LEAD10_TOP * Screen.TwipsPerPixelY
    Lead10.Left = LEAD10_LEFT * Screen.TwipsPerPixelX
    Lead10.Height = LEAD10_HEIGHT * Screen.TwipsPerPixelY
    Lead10.Width = LEAD10_WIDTH * Screen.TwipsPerPixelX
    
    Lead11.Top = LEAD11_TOP * Screen.TwipsPerPixelY
    Lead11.Left = LEAD11_LEFT * Screen.TwipsPerPixelX
    Lead11.Height = LEAD11_HEIGHT * Screen.TwipsPerPixelY
    Lead11.Width = LEAD11_WIDTH * Screen.TwipsPerPixelX
    
    Lead12.Top = LEAD12_TOP * Screen.TwipsPerPixelY
    Lead12.Left = LEAD12_LEFT * Screen.TwipsPerPixelX
    Lead12.Height = LEAD12_HEIGHT * Screen.TwipsPerPixelY
    Lead12.Width = LEAD12_WIDTH * Screen.TwipsPerPixelX
    
    Lead13.Top = LEAD13_TOP * Screen.TwipsPerPixelY
    Lead13.Left = LEAD13_LEFT * Screen.TwipsPerPixelX
    Lead13.Height = LEAD13_HEIGHT * Screen.TwipsPerPixelY
    Lead13.Width = LEAD13_WIDTH * Screen.TwipsPerPixelX
    
    Lead14.Top = LEAD14_TOP * Screen.TwipsPerPixelY
    Lead14.Left = LEAD14_LEFT * Screen.TwipsPerPixelX
    Lead14.Height = LEAD14_HEIGHT * Screen.TwipsPerPixelY
    Lead14.Width = LEAD14_WIDTH * Screen.TwipsPerPixelX
    
    Lead15.Top = LEAD15_TOP * Screen.TwipsPerPixelY
    Lead15.Left = LEAD15_LEFT * Screen.TwipsPerPixelX
    Lead15.Height = LEAD15_HEIGHT * Screen.TwipsPerPixelY
    Lead15.Width = LEAD15_WIDTH * Screen.TwipsPerPixelX
    
    leadProbe.Top = PROBE_TOP * Screen.TwipsPerPixelY
    leadProbe.Left = PROBE_LEFT * Screen.TwipsPerPixelX
    leadProbe.Height = PROBE_HEIGHT * Screen.TwipsPerPixelY
    leadProbe.Width = PROBE_WIDTH * Screen.TwipsPerPixelX
    
    lblProbe.Top = PROBE_TOP * Screen.TwipsPerPixelY        'takeout shape
    lblProbe.Left = PROBE_LEFT * Screen.TwipsPerPixelX      'takeout shape
    lblProbe.Height = PROBE_HEIGHT * Screen.TwipsPerPixelY  'takeout shape
    lblProbe.Width = PROBE_WIDTH * Screen.TwipsPerPixelX    'takeout shape
    
    lblDiff.Top = DIFF_TOP * Screen.TwipsPerPixelY
    lblDiff.Left = DIFF_LEFT * Screen.TwipsPerPixelX
    lblDiff.Height = DIFF_HEIGHT * Screen.TwipsPerPixelY
    lblDiff.Width = DIFF_WIDTH * Screen.TwipsPerPixelX
    
    lblDiffTouch.Top = DIFF_TOP * Screen.TwipsPerPixelY
    lblDiffTouch.Left = (DIFF_LEFT - PROBE_WIDTH + DIFF_WIDTH) * Screen.TwipsPerPixelX
    lblDiffTouch.Height = PROBE_HEIGHT * Screen.TwipsPerPixelY
    lblDiffTouch.Width = PROBE_WIDTH * Screen.TwipsPerPixelX
    
    
    'gsubScreenDimensions Me, 640, 0
    'gsubScreenDimensions Me, 1280, 0
    gsubScreenDimensions Me, 1024, 0

End Sub

Private Sub Form_MouseUp(Button As Integer, Shift As Integer, X As Single, Y As Single)
    frmControl.DetectedFormPeck
End Sub

Private Sub Lead1_MouseDown(Button As Integer, Shift As Integer, X As Single, Y As Single)
    frmControl.DetectedPos2
End Sub
Private Sub Lead1_MouseUp(Button As Integer, Shift As Integer, X As Single, Y As Single)
'    frmControl.DetectedPos2
End Sub
Private Sub Lead2_MouseDown(Button As Integer, Shift As Integer, X As Single, Y As Single)
    frmControl.DetectedPos3
End Sub
Private Sub Lead2_MouseUp(Button As Integer, Shift As Integer, X As Single, Y As Single)
 '   frmControl.DetectedPos3
End Sub
Private Sub Lead3_MouseDown(Button As Integer, Shift As Integer, X As Single, Y As Single)
    frmControl.DetectedPos4
End Sub
Private Sub Lead3_MouseUp(Button As Integer, Shift As Integer, X As Single, Y As Single)
  '  frmControl.DetectedPos4
End Sub
Private Sub Lead4_MouseDown(Button As Integer, Shift As Integer, X As Single, Y As Single)
    frmControl.DetectedPos5
End Sub
Private Sub Lead4_MouseUp(Button As Integer, Shift As Integer, X As Single, Y As Single)
   ' frmControl.DetectedPos5
End Sub
Private Sub Lead5_MouseDown(Button As Integer, Shift As Integer, X As Single, Y As Single)
    frmControl.DetectedPos6
End Sub
Private Sub Lead5_MouseUp(Button As Integer, Shift As Integer, X As Single, Y As Single)
    'frmControl.DetectedPos6
End Sub
Private Sub Lead6_MouseDown(Button As Integer, Shift As Integer, X As Single, Y As Single)
    frmControl.DetectedPos7
End Sub
Private Sub Lead6_MouseUp(Button As Integer, Shift As Integer, X As Single, Y As Single)
    'frmControl.DetectedPos7
End Sub
Private Sub Lead7_MouseDown(Button As Integer, Shift As Integer, X As Single, Y As Single)
    frmControl.DetectedPos8
End Sub
Private Sub Lead7_MouseUp(Button As Integer, Shift As Integer, X As Single, Y As Single)
    'frmControl.DetectedPos8
End Sub
Private Sub Lead8_MouseDown(Button As Integer, Shift As Integer, X As Single, Y As Single)
    frmControl.DetectedPos9
End Sub
Private Sub Lead8_MouseUp(Button As Integer, Shift As Integer, X As Single, Y As Single)
    'frmControl.DetectedPos9
End Sub
Private Sub Lead9_MouseDown(Button As Integer, Shift As Integer, X As Single, Y As Single)
    frmControl.DetectedPos10
End Sub
Private Sub Lead9_MouseUp(Button As Integer, Shift As Integer, X As Single, Y As Single)
    'frmControl.DetectedPos10
End Sub
Private Sub Lead10_MouseDown(Button As Integer, Shift As Integer, X As Single, Y As Single)
    frmControl.DetectedPos11
End Sub
Private Sub Lead10_MouseUp(Button As Integer, Shift As Integer, X As Single, Y As Single)
   ' frmControl.DetectedPos11
End Sub
Private Sub Lead11_MouseDown(Button As Integer, Shift As Integer, X As Single, Y As Single)
    frmControl.DetectedPos12
End Sub
Private Sub Lead11_MouseUp(Button As Integer, Shift As Integer, X As Single, Y As Single)
   ' frmControl.DetectedPos12
End Sub
Private Sub Lead12_MouseDown(Button As Integer, Shift As Integer, X As Single, Y As Single)
    frmControl.DetectedPos13
End Sub
Private Sub Lead12_MouseUp(Button As Integer, Shift As Integer, X As Single, Y As Single)
   ' frmControl.DetectedPos13
End Sub
Private Sub Lead13_MouseDown(Button As Integer, Shift As Integer, X As Single, Y As Single)
    frmControl.DetectedPos14
End Sub
Private Sub Lead13_MouseUp(Button As Integer, Shift As Integer, X As Single, Y As Single)
   ' frmControl.DetectedPos14
End Sub
Private Sub Lead14_MouseDown(Button As Integer, Shift As Integer, X As Single, Y As Single)
    frmControl.DetectedPos15
End Sub
Private Sub Lead14_MouseUp(Button As Integer, Shift As Integer, X As Single, Y As Single)
   ' frmControl.DetectedPos15
End Sub
Private Sub Lead15_MouseDown(Button As Integer, Shift As Integer, X As Single, Y As Single)
    frmControl.DetectedPos16
End Sub
Private Sub Lead15_MouseUp(Button As Integer, Shift As Integer, X As Single, Y As Single)
   ' frmControl.DetectedPos16
End Sub


Private Sub leadProbe_MouseDown(Button As Integer, Shift As Integer, X As Single, Y As Single)
    frmControl.handleReadySignalTouch
End Sub

Private Sub leadStimuli_MouseDown(Button As Integer, Shift As Integer, X As Single, Y As Single)
    frmControl.DetectedPos1
End Sub
Private Sub leadStimuli_MouseUp(Button As Integer, Shift As Integer, X As Single, Y As Single)
   ' frmControl.DetectedPos1
End Sub
