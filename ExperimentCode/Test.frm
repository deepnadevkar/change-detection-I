VERSION 5.00
Object = "{02B97C13-ED1E-11CD-A08B-00AA00575482}#1.0#0"; "mhiinp32.ocx"
Object = "{4DE9E2A3-150F-11CF-8FBF-444553540000}#4.0#0"; "DlxOCX32.ocx"
Begin VB.Form frmTest 
   BackColor       =   &H00404040&
   BorderStyle     =   1  'Fixed Single
   Caption         =   "Test"
   ClientHeight    =   3270
   ClientLeft      =   14040
   ClientTop       =   2430
   ClientWidth     =   3000
   ControlBox      =   0   'False
   LinkTopic       =   "Form1"
   MaxButton       =   0   'False
   MinButton       =   0   'False
   PaletteMode     =   1  'UseZOrder
   ScaleHeight     =   3270
   ScaleWidth      =   3000
   Begin MhinintLibCtl.MhIntInput intPelletTime 
      DataField       =   "PelletTime"
      DataSource      =   "datDefaults"
      Height          =   315
      Left            =   900
      TabIndex        =   10
      Top             =   930
      Width           =   585
      _Version        =   65536
      _ExtentX        =   1032
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Text            =   "3000"
      Max             =   30000
      Min             =   0
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput intJuiceTime 
      DataField       =   "JuiceTime"
      DataSource      =   "datDefaults"
      Height          =   315
      Left            =   900
      TabIndex        =   8
      Top             =   240
      Width           =   585
      _Version        =   65536
      _ExtentX        =   1032
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Text            =   "3000"
      Max             =   30000
      Min             =   0
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin DlsrLib.DriverLINXSR DriverLINXSR1 
      Left            =   2040
      Top             =   240
      _Version        =   262144
      _ExtentX        =   741
      _ExtentY        =   741
      _StockProps     =   64
   End
   Begin VB.Timer tmrReadTouch 
      Enabled         =   0   'False
      Interval        =   100
      Left            =   1650
      Top             =   1380
   End
   Begin VB.Timer tmrInit 
      Enabled         =   0   'False
      Left            =   0
      Top             =   1350
   End
   Begin VB.Timer tmrMetraByte 
      Enabled         =   0   'False
      Left            =   780
      Top             =   1410
   End
   Begin VB.CommandButton cmdPellets 
      Caption         =   "Hopper"
      Height          =   525
      Left            =   60
      TabIndex        =   1
      Top             =   840
      Width           =   825
   End
   Begin VB.CommandButton cmdJuice 
      Caption         =   "House Light"
      Height          =   525
      Left            =   60
      TabIndex        =   0
      Top             =   150
      Width           =   825
   End
   Begin VB.Label Label34 
      BackColor       =   &H00404040&
      Caption         =   "msec."
      ForeColor       =   &H00E0E0E0&
      Height          =   255
      Left            =   1530
      TabIndex        =   11
      Top             =   990
      Width           =   495
   End
   Begin VB.Label Label35 
      BackColor       =   &H00404040&
      Caption         =   "msec."
      ForeColor       =   &H00E0E0E0&
      Height          =   255
      Left            =   1530
      TabIndex        =   9
      Top             =   300
      Width           =   495
   End
   Begin VB.Label lblMsg 
      Alignment       =   2  'Center
      BorderStyle     =   1  'Fixed Single
      Height          =   285
      Left            =   210
      TabIndex        =   7
      Top             =   2040
      Width           =   1665
   End
   Begin VB.Label Label5 
      Alignment       =   2  'Center
      Caption         =   "Y"
      BeginProperty Font 
         Name            =   "MS Sans Serif"
         Size            =   9.75
         Charset         =   0
         Weight          =   700
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      Height          =   255
      Left            =   1320
      TabIndex        =   6
      Top             =   2340
      Visible         =   0   'False
      Width           =   285
   End
   Begin VB.Label Label4 
      Alignment       =   2  'Center
      Caption         =   "X"
      BeginProperty Font 
         Name            =   "MS Sans Serif"
         Size            =   9.75
         Charset         =   0
         Weight          =   700
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      Height          =   255
      Left            =   480
      TabIndex        =   5
      Top             =   2340
      Visible         =   0   'False
      Width           =   285
   End
   Begin VB.Label lblY 
      Alignment       =   2  'Center
      BorderStyle     =   1  'Fixed Single
      Caption         =   "000"
      BeginProperty Font 
         Name            =   "MS Sans Serif"
         Size            =   13.5
         Charset         =   0
         Weight          =   700
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      Height          =   435
      Left            =   1050
      TabIndex        =   4
      Top             =   2610
      Visible         =   0   'False
      Width           =   825
   End
   Begin VB.Label lblX 
      Alignment       =   2  'Center
      BorderStyle     =   1  'Fixed Single
      Caption         =   "000"
      BeginProperty Font 
         Name            =   "MS Sans Serif"
         Size            =   13.5
         Charset         =   0
         Weight          =   700
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      Height          =   435
      Left            =   210
      TabIndex        =   3
      Top             =   2610
      Visible         =   0   'False
      Width           =   825
   End
   Begin VB.Label Label1 
      Alignment       =   2  'Center
      BorderStyle     =   1  'Fixed Single
      Caption         =   "Inoperative"
      Height          =   1635
      Left            =   60
      TabIndex        =   2
      Top             =   1560
      Width           =   1935
   End
   Begin VB.Menu mnuFIle 
      Caption         =   "&File"
      Begin VB.Menu mnuClose 
         Caption         =   "&Close"
      End
   End
End
Attribute VB_Name = "frmTest"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = False
Option Explicit

'Const LogicalDevice As Integer = 0  'DriverLINX Logical Device
'Private LogicalChannel As Integer   'DriverLINX Logical Channel
'Const gain As Single = 0            'DriverLINX Gain
'Const BackGroundForeGround As Integer = 0     'IRQ or DMA = 1, Polled = 0
'Const ChannelGain As Single = -1#   'Gain setting for Logical Channel

Dim xmouse As Single
Dim ymouse As Single

Private Sub cmdJuice_Click()
   ' old code
   ' agOutp METRABYTE_PORT_B, APPARATUS_HOUSE_LIGHT
    ' setup DO Service Request for polled mode Output operation

With DriverLINXSR1
    .Req_subsystem = DL_DO
    .Req_mode = DL_POLLED
    .Req_op = DL_START
    .Evt_Str_type = DL_NULLEVENT
    .Evt_Stp_type = DL_NULLEVENT
    .Evt_Tim_type = DL_NULLEVENT
    .Sel_chan_format = DL_tNATIVE
    .Sel_chan_N = 1
    .Sel_chan_start = 1   'Channel 0 is configured for Output already
    .Sel_chan_startGainCode = 0
    .Sel_buf_samples = 1
    .Sel_buf_N = 0
    .Refresh
    .Res_Sta_ioValue = APPARATUS_HOUSE_LIGHT
    .Refresh
    Label1.Caption = .Message
End With

    tmrMetraByte.Interval = intJuiceTime
    tmrMetraByte.Enabled = True

End Sub

Private Sub cmdPellets_Click()

    'agOutp METRABYTE_PORT_B, APPARATUS_HOPPER
    'agOutp METRABYTE_PORT_B, APPARATUS_HOPPER_LIGHT
    'agOutp METRABYTE_PORT_B, APPARATUS_HOPPER + APPARATUS_HOPPER_LIGHT
    
    With DriverLINXSR1
    .Req_subsystem = DL_DO
    .Req_mode = DL_POLLED
    .Req_op = DL_START
    .Evt_Str_type = DL_NULLEVENT
    .Evt_Stp_type = DL_NULLEVENT
    .Evt_Tim_type = DL_NULLEVENT
    .Sel_chan_format = DL_tNATIVE
    .Sel_chan_N = 1
    .Sel_chan_start = 1   'Channel 0 is configured for Output already
    .Sel_chan_startGainCode = 0
    .Sel_buf_samples = 1
    .Sel_buf_N = 0
    .Refresh
    .Res_Sta_ioValue = APPARATUS_PELLETS
    .Refresh
    Label1.Caption = .Message
    End With
    
    
    tmrMetraByte.Interval = intPelletTime
    tmrMetraByte.Enabled = True

End Sub


Private Sub Command1_Click()

Dim dummy As Integer
dummy = 1
Open "\\box4\c\jeff\jeff.txt" For Output As #1
        For dummy = 1 To 5
           Print #1, dummy; " JSK "
            'Debug.Print i,
            
        Next
            
    Close #1
End Sub

Private Sub Form_Activate()

    Refresh
    
    
     
    'lblMsg.Caption = "Initializing"
    'gsInitState = "Init"
    'tmrInit.Interval = 100
    'tmrInit.Enabled = True
    
End Sub

Private Sub Form_Load()
    
    'initialize PCI driver
    With DriverLINXSR1
        .Req_device = 0
        .Req_mode = DL_OTHER
        .Req_op = DL_INITIALIZE
        .Req_DLL_name = "KPCIPIO"
        .Refresh
        Label1.Caption = .Message
    End With
    
    ' code below will Configure Channel 1 for Output
    With DriverLINXSR1
        .Req_subsystem = DL_DO  ' the subsystem to assign the channel to
        .Req_mode = DL_OTHER
        .Req_op = DL_CONFIGURE  ' it is a configuration type operation
        .Evt_Tim_type = DL_DIOSETUP
        .Evt_Tim_dioChannel = 1   ' configure channel 1 or Port b as Output
        .Evt_Tim_dioMode = DL_DIO_BASIC
        .Evt_Str_type = DL_NULLEVENT
        .Evt_Stp_type = DL_NULLEVENT
        .Sel_chan_N = 0
        .Refresh
        lblMsg.Caption = .Message
    End With
    
    
'old metrabyte
    'agOutp METRABYTE_CONTROL, APPARATUS_SET_B_OUTPUT

End Sub


Private Sub Form_MouseMove(Button As Integer, Shift As Integer, X As Single, Y As Single)

xmouse = X
ymouse = Y
End Sub

Private Sub mnuClose_Click()

   ' gsubDisableTouchscreen tmrReadTouch
    DriverLINXSR1.Req_DLL_name = "" 'unloads driver
    Unload Me
    
End Sub


Private Sub tmrInit_Timer()
    
  '  gsubInitTouchscreen tmrInit, gsInitState, tmrReadTouch
    tmrReadTouch.Interval = 100
    tmrReadTouch.Enabled = True

End Sub

Private Sub tmrMetraByte_Timer()
    
    tmrMetraByte.Enabled = False
    tmrMetraByte.Interval = 0
    
   'agOutp METRABYTE_PORT_B, APPARATUS_OFF
    DriverLINXSR1.Res_Sta_ioValue = APPARATUS_OFF
    DriverLINXSR1.Refresh
    
    
End Sub


Private Sub tmrReadTouch_Timer()

    Dim iX As Integer
    Dim iY As Integer
    
    lblMsg.Caption = "Ready"

    'gsubGetXY iX, iY
    lblX.Caption = xmouse
    lblY.Caption = ymouse
    'lblX.Caption = iX
    'lblY.Caption = iY

End Sub


