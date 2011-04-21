#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <stdlib.h>
#include <iostream>
using namespace std;

#include <QMainWindow>
#include <QFileDialog>
#include <QtSql>

MainWindow::MainWindow(QWidget *parent) :    QMainWindow(parent),ui(new Ui::MainWindow){
	ui->setupUi(this);
}

MainWindow::~MainWindow(){
	delete ui;
}

/*
  could allow for multiple connections but not doing that for now
 */
void MainWindow::load_sequencesqlite(){
	if (sdb.connectionName().length() != 0){
		cout << "existing connection" << endl;
		sdb.close();
		sdb.removeDatabase(sdb.connectionName());
	}
	QString filename = QFileDialog::getOpenFileName( this, tr("Open Document"), QDir::currentPath(), tr("Database files (*.db *.sql *.sqlite);;All files (*.*)"), 0, QFileDialog::DontUseNativeDialog );
	cout << filename.toStdString() << endl;
	sdb = QSqlDatabase::addDatabase("QSQLITE");
	sdb.setDatabaseName(filename);
	if (!sdb.open())
			cout << "open failed" << endl;
	else
		ui->statusBar->showMessage("sequence database opened");

//	QSqlQuery q;
//	q.exec("SELECT * from taxonomy where edited_name = 'Viburnum'");
//	if (q.next()){
//		cout << q.value(0).toString().toStdString() << endl;
//	}

}

void MainWindow::get_sequencedb_stats(){

}
