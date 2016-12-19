unit Unit1;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, Forms, Controls, Graphics, Dialogs, StdCtrls,
  ExtCtrls, Math;

type

  { TForm1 }

  TForm1 = class(TForm)
    CheckBox1: TCheckBox;
    czyscBut: TButton;
    bEdit: TEdit;
    bLabel: TLabel;
    Edit1: TEdit;
    symulacjaOptBut: TButton;
    OptymalizacjaCheckBox: TCheckBox;
    optymalizacjaPlikEdit: TEdit;
    Label6: TLabel;
    symulacjaOptWykres: TShape;
    xminImpLabel: TLabel;
    xmaxImpLabel: TLabel;
    vxImpEdit: TEdit;
    vxImpLabel: TLabel;
    vminImpEdit: TEdit;
    vmaxImpEdit: TEdit;
    vminImpLabel: TLabel;
    vmaxImpLabel: TLabel;
    LiniowaButton: TButton;
    OptymalizujBut: TButton;
    ImpulsCheckBox: TCheckBox;
    ImpulsPlikEdit: TEdit;
    impulsBut: TButton;
    dtEdit: TEdit;
    dtLabel: TLabel;
    impulsWykres: TShape;
    Label1: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    wymuszenieWykres: TShape;
    WskaznikListBox: TListBox;
    vqtrLabel: TLabel;
    xqtrLabel: TLabel;
    vq_trLabel1: TLabel;
    vuLabel2: TLabel;
    trLabel2: TLabel;
    xq_trLabel: TLabel;
    vq_trLabel: TLabel;
    xq_trLabel1: TLabel;
    xvEdit: TEdit;
    xvLabel: TLabel;
    vminEdit: TEdit;
    vmaxEdit: TEdit;
    vminLabel: TLabel;
    vmaxLabel: TLabel;
    xminLabel: TLabel;
    xmaxLabel: TLabel;
    StopBut: TButton;
    symulacjaBut: TButton;
    cqLabel: TLabel;
    tsymEdit: TEdit;
    tsymLabel: TLabel;
    symulacjaWykres: TShape;
    vuEdit: TEdit;
    vuLabel: TLabel;
    tuEdit: TEdit;
    tuLabel: TLabel;
    trEdit: TEdit;
    fQEdit: TEdit;
    fQLabel: TLabel;
    trLabel: TLabel;
    LEdit: TEdit;
    mQEdit: TEdit;
    mQLabel: TLabel;
    lLabel: TLabel;
    procedure czyscButClick(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure impulsButClick(Sender: TObject);
    procedure LiniowaButtonClick(Sender: TObject);
    procedure OptymalizujButClick(Sender: TObject);
    procedure StopButClick(Sender: TObject);
    procedure symulacjaButClick(Sender: TObject);
    procedure symulacjaOptButClick(Sender: TObject);
    procedure wymuszenie();
    procedure rysujWykresy();
    procedure rysujImpuls();
    procedure rysujWymuszenie();
    procedure rysujSymulacjaOpt();

  private
    { private declarations }
  public
    { public declarations }
  end;

var
  Form1: TForm1;
  mq : real = 1000; // masa ladunku [kg]
  L : real = 5; // dlugosc lin [m]
  cq : real; // sztywnosc pozioma zawieszonego ladunku [N/m]
  fq : real = 0; // opory powietrza
  tr : real = 3; // czas rozruchu [s]
  tu : real = 10; // czas jazdy ustalonej [s]
  vu : real = 10; // predkosc ustalona jazdy [m/min]
  tsym : real = 30; // czas symulacji
  t, v, th : real;
  dt : real = 0.001;
  vq, dvq, vs, xq, dxq : real;
  vq_tr, xq_tr : real;
  vqtr, xqtr : real;

  // zmienne do wykresu
  szer : integer = 450;
  wys : integer = 300;
  szer1 : integer = 300;
  wys1 : integer = 180;
  skx, skv, skt, skv1, skt1,sktImp, skxImp, skvImp : real; // skale
  xmax, xmin : real;
  vmax : real = 20; // [m/min]
  vmin : real = -20; // [m/min]
  xv : real = 0.01; // [-] stosunek skali x do skali v
  vxImp : real = 100; // [-] stosunek skali v do x dla Impulsu
  vmaxImp : real = 1; // [m/min]
  vminImp : real = -1;
  xmaxImp, xminImp : real;

  // zmienne do impulsu
  xqTablicaImp : Array[0..10000] of Real;
  vqTablicaImp : Array[0..10000] of Real;
  vsTablica : Array[0..10000] of Real;
  i, j, k : LongInt;

  // optymalizacja
  b : real = 0.2; // wsp. beta
  Iq :real; // wsp. jakosci
  xqu, vqu : real;
  p : Array[0..10000] of Real;
  vsTablica2 : Array[0..10000] of Real;
  amax: real =0.167; // m/s^2
  deltav : real;
  znak : Integer;

  // zapis do pliku
  plik:text; nazwa: string;
  licznik:Integer;

const
  g = 9.81;

implementation

{$R *.lfm}

{ TForm1 }

procedure otw_plik;
  begin
    assignfile (plik,nazwa);
    rewrite(plik);
    //writeln(plik,'t;vs;xq;vq;i');
  end;

procedure TForm1.FormCreate(Sender: TObject);
begin
  mQLabel.Caption:='mQ = '+floattostr(mq);
  mQEdit.Text:=floattostr(mq);
  lLabel.Caption:='l = '+floattostr(L);
  lEdit.Text:=floattostr(L);
  cQLabel.Caption:='cQ = '+floattostr(cq);
  fQLabel.Caption:='fQ = '+floattostr(fq);
  fQEdit.Text:=floattostr(fq);
  trLabel.Caption:='tr = '+floattostr(tr);
  trEdit.Text:=floattostr(tr);
  tuLabel.Caption:='tu = '+floattostr(tu);
  tuEdit.Text:=floattostr(tu);
  vuLabel.Caption:='vu = '+floattostr(vu);
  vuEdit.Text:=floattostr(vu);
  tsymLabel.Caption:='tsym = '+floattostr(tsym);
  tsymEdit.Text:=floattostr(tsym);
  dtLabel.Caption:='dt = '+floattostr(dt);
  dtEdit.Text:=floattostr(dt);
  xvEdit.Text:=floattostr(xv);
  vmaxEdit.Text:=floattostr(vmax);
  vminEdit.Text:=floattostr(vmin);
  vxImpEdit.Text:=floattostr(vxImp);
  vmaxImpEdit.Text:=floattostr(vmaxImp);
  vminImpEdit.Text:=floattostr(vminImp);
  bEdit.Text:=floattostr(b);
  // wykres
  symulacjaWykres.Height:=wys;
  symulacjaWykres.Width:=szer;
  impulsWykres.Height:=wys;
  impulsWykres.Width:=szer;
  wymuszenieWykres.Height:=wys1;
  wymuszenieWykres.Width:=szer1;
  symulacjaOptWykres.Height:=wys;
  symulacjaOptWykres.Width:=szer;
  // zapis do pliku
  ImpulsPlikEdit.Text:='impuls.txt';
  optymalizacjaPlikEdit.Text:='optymalizacja.txt'
end;

procedure TForm1.czyscButClick(Sender: TObject);
begin
  wymuszenieWykres.Canvas.Rectangle(0,0,szer1,wys1);
  wskaznikListBox.Clear;
end;


procedure TForm1.impulsButClick(Sender: TObject);
begin
  // wczytanie zmiennych
  mq := strtofloat(mqEdit.text);
  mQLabel.Caption:='mQ = '+floattostr(mq);
  L := strtofloat(lEdit.text);
  lLabel.Caption:='l = '+floattostr(l);
  fq := strtofloat(fqEdit.text);
  fQLabel.Caption:='fQ = '+floattostr(fq);
  tr := strtofloat(trEdit.text);
  trLabel.Caption:='tr = '+floattostr(tr);
  tu := strtofloat(tuEdit.text);
  tuLabel.Caption:='tu = '+floattostr(tu);
  vu := strtofloat(vuEdit.text);
  vuLabel.Caption:='vu = '+floattostr(vu);
  tsym := strtofloat(tsymEdit.text);
  tsymLabel.Caption:='tsym = '+floattostr(tsym);
  dt := strtofloat(dtEdit.text);
  dtLabel.Caption:='dt = '+floattostr(dt);
  vxImp := strtofloat(vxImpEdit.text);
  vmaxImp := strtofloat(vmaxImpEdit.text);
  vminImp := strtofloat(vminImpEdit.text);

  // obliczenia
  skvImp:=(wys)/(vmaxImp+abs(vminImp));
  sktImp:=(szer)/(tr);
  xmaxImp:=vmaxImp/vxImp;
  xminImp:=vminImp/vxImp;
  skxImp:=(wys)/(xmaxImp+abs(xminImp));
  xmaxImpLabel.Caption:='xmax='+floattostr(xmaxImp);
  xminImpLabel.Caption:='xmin='+floattostr(xminImp);
  cq:=mq*g/L;
  cQLabel.Caption:='cQ = '+floattostr(cq);

  // czyszczenie wykresu
  impulsWykres.Canvas.Rectangle(0,0,szer,wys);
  impulsWykres.Canvas.MoveTo(0,round(vmaxImp*skvImp));
  impulsWykres.Canvas.LineTo(szer,round(vmaxImp*skvImp));

  // zapis do pliku
  if impulscheckbox.Checked then
    begin
      nazwa:=ImpulsPlikEdit.text;
      otw_plik;
      writeln(plik,'t;vs;xq;vq;i');
    end;

  // warunki początkowe
  xq:= 0; vq:= 0; t:= 0; i:=0;

  // pętla symulacji
  repeat
    // impuls
    if ((t >= 0) and (t < dt) ) then v := 1 else v := 0;
    // symulacja
    vs:=v;
    dvq:=-(fq/mq)*vq+(cq/mq)*xq;
    dxq:=vs-vq;
    xqTablicaImp[i]:=xq;
    vqTablicaImp[i]:=vq;
    rysujImpuls();
    if impulscheckbox.Checked then
      writeln(plik,t:6:3,';',vs:6:3,';',xq:9:6,';',vq:9:6,';',i:5);
    vq:=vq+dvq*dt;
    xq:=xq+dxq*dt;
    t:=t+dt;
    i:=i+1;
  until t>(tr+dt/2);
  // zapis do pliku
  if impulscheckbox.checked then closefile(plik);

end;

procedure TForm1.LiniowaButtonClick(Sender: TObject);
begin
  // oblicza skale do wykresu
  skt1:= (szer1)/(tr);
  skv1:= (wys1)/(vu/60);
  // stwarza wektor sterowania początkowego - funkcja liniowa
  t:=0;
  for j:=0 to i do
      begin
        vs:=((vu/60)/tr)*t;
        vsTablica[j]:=vs;
        wymuszenieWykres.Canvas.Pixels[round(t*skt1), round(skv1*(vu/60-vsTablica[j]))]:=clgreen;
        t:=t+dt;
      end;
  trLabel2.Caption:='tr = '+floattostr(t);
  vuLabel2.Caption:='vu = '+floattostr(vsTablica[i]*60);
  // obliczenie wartosci zmiennych stanu po skonczonym sterowaniu - funkcja liniowa
  k:=round(tr*1000); xqtr:=0; vqtr:=0; t:=0;
  for j:=0 to i do
      begin
        xqtr:= xqtr+(xqTablicaImp[k-j]*vsTablica[j]);
        vqtr:= vqtr+(vqTablicaImp[k-j]*vsTablica[j]);
        t:=t+dt;
      end;
  xqtrLabel.Caption:=floattostr(xqtr);
  vqtrLabel.Caption:=floattostr(vqtr);
  // wskaźnik jakosci
  Iq := 0.5*mq*(vqtr-vqu)*(vqtr-vqu)+0.5*cq*(xqtr-xqu)*(xqtr-xqu);
  wskazniklistbox.Items.Add('Iq='+floattostr(Iq));
  // rysuje wykres
  {t := 0;
  for j:=0 to i do
  begin
  wymuszenieWykres.Canvas.Pixels[round(t*skt1), round(skv1*(vu/60-vsTablica[j]))]:=clgreen;
  t:=t+dt;}
  end;


//end;

procedure TForm1.OptymalizujButClick(Sender: TObject);

var
  IqOPT : real;
  exit : boolean;
  alfa : real;
begin
  // wczytanie b
  b:=strtofloat(bEdit.text);
  bLabel.Caption:='beta='+floattostr(b);
  // obliczenie wektora p
  k:=round(tr*1000);
  deltav:= amax*dt;
  exit:=False;
  IqOPT:=Iq;
  //glowna petla
  repeat
    bLabel.Caption:='beta='+floattostr(b);
  for j:=0 to i do
        p[j]:=(-b)*((mq*(vqtr-vqu)*vqTablicaImp[k-j])+(cq*(xqtr-xqu)*xqTablicaImp[k-j]));
  // stworzenie wektora sterowania vs
  {if optymalizacjacheckbox.Checked then
    begin
      nazwa:=optymalizacjaPlikEdit.text;
      otw_plik;
      writeln(plik,'t;vs');
    end; }

  t:=0; //k:=round(tr*1000);
  for j:=0 to i do
      begin
        vsTablica2[j]:=vsTablica[j]+p[j];
        if ((j>0) and (j<i)) then
          begin
            if (abs(vsTablica2[j]-vsTablica2[j-1])> deltav) then
              begin
                znak:= sign(vsTablica2[j]-vsTablica2[j-1]);
                vsTablica2[j]:=vsTablica2[j-1]+deltav*znak;
              end;
            if ((((vu/60)-vsTablica2[j])/(tr-t))>amax) then
              vsTablica2[j]:=vsTablica2[j-1]+deltav;
            if vsTablica2[j]>(vu/60)then vsTablica2[j]:=vu/60;
            if vsTablica2[j]<0 then vsTablica2[j]:=0;
          end;
        vsTablica2[0]:=0;
        vsTablica2[i]:=vu/60;
        //vsTablica[j]:=vsTablica2[j];
        //if optymalizacjacheckbox.Checked then
           //writeln(plik,t:6:3,';',vsTablica[j]:6:3,';',j:5);
        rysujWymuszenie();
        t:=t+dt;
      end;
  {k:=3000;
  for j:=0 to k do
      begin
        if (j<3000) then
            if (abs(vsTablica2[k-j]-vsTablica2[k-j-1])> deltav) then
              vsTablica2[k-j-1]:=vsTablica2[k-j]-deltav;
      end;
  t:=0; k:=3000;
  for j:=0 to k do
      begin
      vsTablica[j]:=vsTablica2[j];
  if optymalizacjacheckbox.Checked then
           writeln(plik,t:6:3,';',vsTablica[j]:6:3,';',j:5);
        rysujWymuszenie();
        t:=t+dt;
      end; }
  // zapis do pliku
  //if optymalizacjacheckbox.checked then closefile(plik);
  // obliczenie wartosci zmiennych stanu po skonczonym sterowaniu
  k:=round(tr*1000); xqtr:=0; vqtr:=0; t:=0;
  for j:=0 to i do
      begin
        xqtr:= xqtr+(xqTablicaImp[k-j]*vsTablica2[j]);
        vqtr:= vqtr+(vqTablicaImp[k-j]*vsTablica2[j]);
        t:=t+dt;
      end;
  xqtrLabel.Caption:=floattostr(xqtr);
  vqtrLabel.Caption:=floattostr(vqtr);
  // wskaźnik jakosci
  Iq := 0.5*mq*(vqtr-vqu)*(vqtr-vqu)+0.5*cq*(xqtr-xqu)*(xqtr-xqu);
  wskazniklistbox.Items.Add('Iq='+floattostr(Iq));
  if Iq>IqOPT then b:=b/2;
  if Iq<=IqOPT then
    begin
      alfa:=IqOPT-Iq;
      if ((alfa)<power(10,-30)) then exit:=true;
      for j:=0 to i do
      vsTablica[j]:=vsTablica2[j];
      IqOPT:=Iq;
    end;
  until (exit=true);
  wymuszenieWykres.Canvas.Rectangle(0,0,szer1,wys1);
  t := 0;
  for j:=0 to i do
  begin
  wymuszenieWykres.Canvas.Pixels[round(t*skt1), round(skv1*(vu/60-vsTablica[j]))]:=clgreen;
  t:=t+dt;
  end;
end;


procedure TForm1.wymuszenie();
begin
  if (t<=1) then v := 0;
  if ((t>1) and (t<=(tr+1))) then v := vu/tr * (t-1);
  if ((t>(tr+1)) and (t<=(tu+tr+1))) then v := vu;
  if ((t>(tu+tr+1)) and (t<=(th+tu+tr+1))) then v := vu - vu/th * (t-(tu+tr+1));
  if (t>(th+tu+tr+1)) then v:=0;
end;

procedure TForm1.rysujWykresy();
begin
  symulacjaWykres.Canvas.Pixels[round(t*skt), round(skv*(vmax-v))]:=clblue;
  symulacjaWykres.Canvas.Pixels[round(t*skt), round(skv*(vmax-vq*60))]:=clgreen;
  symulacjaWykres.Canvas.Pixels[round(t*skt), round(skx*(xmax-xq))]:=clred;
end;

procedure TForm1.rysujImpuls();
begin
  impulsWykres.Canvas.Pixels[round(t*sktImp), round(skvImp*(vmaxImp-v*60))]:=clblue;
  impulsWykres.Canvas.Pixels[round(t*sktImp), round(skvImp*(vmaxImp-vq*60))]:=clgreen;
  impulsWykres.Canvas.Pixels[round(t*sktImp), round(skxImp*(xmaxImp-xq))]:=clred;
end;

procedure TForm1.rysujWymuszenie();
begin
  wymuszenieWykres.Canvas.Pixels[round(t*skt1), round(skv1*(vu/60-vsTablica2[j]))]:=clgreen;
end;

procedure TForm1.rysujSymulacjaOpt();
begin
  symulacjaOptWykres.Canvas.Pixels[round(t*skt), round(skv*(vmax-vs*60))]:=clblue;
  symulacjaOptWykres.Canvas.Pixels[round(t*skt), round(skv*(vmax-vq*60))]:=clgreen;
  symulacjaOptWykres.Canvas.Pixels[round(t*skt), round(skx*(xmax-xq))]:=clred;
end;

procedure TForm1.StopButClick(Sender: TObject);
begin
  halt;
end;

procedure TForm1.symulacjaButClick(Sender: TObject);
begin

  // wczytanie zmiennych
  mq := strtofloat(mqEdit.text);
  mQLabel.Caption:='mQ = '+floattostr(mq);
  L := strtofloat(lEdit.text);
  lLabel.Caption:='l = '+floattostr(l);
  fq := strtofloat(fqEdit.text);
  fQLabel.Caption:='fQ = '+floattostr(fq);
  tr := strtofloat(trEdit.text);
  trLabel.Caption:='tr = '+floattostr(tr);
  tu := strtofloat(tuEdit.text);
  tuLabel.Caption:='tu = '+floattostr(tu);
  vu := strtofloat(vuEdit.text);
  vuLabel.Caption:='vu = '+floattostr(vu);
  tsym := strtofloat(tsymEdit.text);
  tsymLabel.Caption:='tsym = '+floattostr(tsym);
  dt := strtofloat(dtEdit.text);
  dtLabel.Caption:='dt = '+floattostr(dt);
  xv := strtofloat(xvEdit.text);
  vmax := strtofloat(vmaxEdit.text);
  vmin := strtofloat(vminEdit.text);

  // obliczenia
  skv:=(wys)/(vmax+abs(vmin));
  skt:=(szer)/tsym;
  xmax:=vmax*xv;
  xmin:=vmin*xv;
  skx:=(wys)/(xmax+abs(xmin));
  xmaxLabel.Caption:='xmax='+floattostr(xmax);
  xminLabel.Caption:='xmin='+floattostr(xmin);
  cq:=mq*g/L;
  cQLabel.Caption:='cQ = '+floattostr(cq);
  // warunki ustalone
  vqu := vu/60;
  xqu := (fq/cq)*vqu;

  // czyszczenie wykresu
  symulacjaWykres.Canvas.Rectangle(0,0,szer,wys);
  symulacjaWykres.Canvas.MoveTo(0,round(vmax*skv));
  symulacjaWykres.Canvas.LineTo(szer,round(vmax*skv));

  // warunki początkowe
  vq:=0; xq:=0;
  th:=tr;
  t:=0;
  // pętla symulacji
  repeat
    wymuszenie();
    vs:=v/60;
    dvq:=-(fq/mq)*vq+(cq/mq)*xq;
    dxq:=vs-vq;
    rysujWykresy();
    if ((t>=(tr+1)) and (t<(tr+1+dt))) then
       begin
         xq_tr:= xq;
         vq_tr:= vq;
       end;
    vq:=vq+dvq*dt;
    xq:=xq+dxq*dt;
    t:=t+dt;
  until t>=tsym;
  xq_trLabel.Caption:=floattostr(xq_tr);
  vq_trLabel.Caption:=floattostr(vq_tr);

end;

procedure TForm1.symulacjaOptButClick(Sender: TObject);
begin
  // czyszczenie wykresu
  symulacjaOptWykres.Canvas.Rectangle(0,0,szer,wys);
  symulacjaOptWykres.Canvas.MoveTo(0,round(vmax*skv));
  symulacjaOptWykres.Canvas.LineTo(szer,round(vmax*skv));

  // zapis do pliku
  if checkbox1.Checked then
    begin
      nazwa:=Edit1.text;
      otw_plik;
      writeln(plik,'t;vs;xq;vq');
    end;
  // warunki początkowe
  vq:= 0; xq:=0; t:=0; j:=0;
  // pętla symulacji
  repeat
    if ((j>=0) and (j<=i)) then vs:= vsTablica[j] else vs:= vu/60;
    dvq:=-(fq/mq)*vq+(cq/mq)*xq;
    dxq:=vs-vq;
    if checkbox1.Checked then
      writeln(plik,t:6:3,';',vs:6:3,';',xq:9:6,';',vq);
    rysujSymulacjaOpt();
    vq:=vq+dvq*dt;
    xq:=xq+dxq*dt;
    t:=t+dt;
    j:=j+1;
  until t>=tsym ;
  // zapis do pliku
  if checkbox1.checked then closefile(plik);
end;

end.

