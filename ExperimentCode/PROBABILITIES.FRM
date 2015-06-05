VERSION 5.00
Object = "{02B97C13-ED1E-11CD-A08B-00AA00575482}#1.0#0"; "mhiinp32.ocx"
Begin VB.Form frmProbabilities 
   BackColor       =   &H00404040&
   BorderStyle     =   1  'Fixed Single
   ClientHeight    =   1905
   ClientLeft      =   15
   ClientTop       =   2025
   ClientWidth     =   2925
   ControlBox      =   0   'False
   LinkTopic       =   "Form1"
   MaxButton       =   0   'False
   MinButton       =   0   'False
   PaletteMode     =   1  'UseZOrder
   ScaleHeight     =   1905
   ScaleWidth      =   2925
   Begin MhinintLibCtl.MhIntInput intPelletProb 
      DataField       =   "PelletProbability"
      DataSource      =   "datDefaults"
      Height          =   315
      Left            =   1740
      TabIndex        =   1
      Top             =   660
      Width           =   555
      _Version        =   65536
      _ExtentX        =   979
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Max             =   100
      Min             =   0
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin MhinintLibCtl.MhIntInput intJuiceProb 
      DataField       =   "JuiceProbability"
      DataSource      =   "datDefaults"
      Height          =   315
      Left            =   1740
      TabIndex        =   3
      Top             =   330
      Visible         =   0   'False
      Width           =   555
      _Version        =   65536
      _ExtentX        =   979
      _ExtentY        =   556
      _StockProps     =   77
      BackColor       =   -2147483633
      BevelSize       =   2
      BorderColor     =   8421504
      BorderStyle     =   0
      FillColor       =   16777215
      LightColor      =   14737632
      ShadowColor     =   4210752
      Max             =   100
      Min             =   0
      AutoHScroll     =   -1  'True
      CaretColor      =   -2147483642
   End
   Begin VB.Data datDefaults 
      Caption         =   "Defaults"
      Connect         =   "Access"
      DatabaseName    =   ""
      DefaultCursorType=   0  'DefaultCursor
      DefaultType     =   2  'UseODBC
      Exclusive       =   0   'False
      Height          =   300
      Left            =   0
      Options         =   0
      ReadOnly        =   0   'False
      RecordsetType   =   1  'Dynaset
      RecordSource    =   ""
      Top             =   0
      Visible         =   0   'False
      Width           =   1965
   End
   Begin VB.CommandButton cmdCancel 
      Caption         =   "Cancel"
      Height          =   345
      Left            =   1560
      TabIndex        =   7
      Top             =   1320
      Width           =   1245
   End
   Begin VB.CommandButton cmdOK 
      Caption         =   "OK"
      Height          =   345
      Left            =   120
      TabIndex        =   6
      Top             =   1320
      Width           =   1245
   End
   Begin VB.Label Label21 
      BackColor       =   &H00404040&
      Caption         =   "Juice Probability:"
      ForeColor       =   &H00FFFFFF&
      Height          =   255
      Left            =   330
      TabIndex        =   5
      Top             =   360
      Visible         =   0   'False
      Width           =   1305
   End
   Begin VB.Label Label22 
      BackColor       =   &H00404040&
      Caption         =   "Pellet Probability:"
      ForeColor       =   &H00FFFFFF&
      Height          =   255
      Left            =   240
      TabIndex        =   4
      Top             =   720
      Width           =   1305
   End
   Begin VB.Label Label26 
      BackColor       =   &H00404040&
      Caption         =   "%"
      ForeColor       =   &H00FFFFFF&
      Height          =   255
      Left            =   2340
      TabIndex        =   2
      Top             =   390
      Visible         =   0   'False
      Width           =   345
   End
   Begin VB.Label Label27 
      BackColor       =   &H00404040&
      Caption         =   "%"
      ForeColor       =   &H00FFFFFF&
      Height          =   255
      Left            =   2340
      TabIndex        =   0
      Top             =   720
      Width           =   345
   End
End
Attribute VB_Name = "frmProbabilities"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = False
Option Explicit

Dim bJuicePelletsChange As Boolean
Private Sub cmdCancel_Click()

    gbRecalculateJuicePellets = False
    
    Unload Me

End Sub

Private Sub cmdOK_Click()
    
    Dim bError As Boolean
    
    bError = False
    
    If Int(intPelletProb.Text) > 100 Or Int(intPelletProb.Text) < 0 Then
        bError = True
        MsgBox "Pellet Probability must between 0 and 100", vbExclamation, "Juice and Pellets Probability"
    End If
    
    If bError = False Then
    
        If bJuicePelletsChange = True Then
        
            datDefaults.Recordset.Edit
            datDefaults.Recordset("JuiceProbability") = 100 - Int(intPelletProb.Text)
            datDefaults.Recordset("PelletProbability") = intPelletProb.Text
            datDefaults.Recordset.Update
            
            gbRecalculateJuicePellets = True
        Else
            gbRecalculateJuicePellets = False
        End If
        
        Unload Me
        
    End If
    
End Sub


Private Sub Form_Activate()
    
    Refresh
    
    datDefaults.DatabaseName = gsProfile
    datDefaults.RecordSource = "SessionInfo"
    datDefaults.Refresh

End Sub

Private Sub Form_Load()

    gbRecalculateJuicePellets = False
    bJuicePelletsChange = False
    
End Sub

Private Sub intJuiceProb_Change()

    bJuicePelletsChange = True
    
End Sub


Private Sub intPelletProb_Change()
    
    bJuicePelletsChange = True

End Sub


