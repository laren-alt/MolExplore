from flask_wtf import FlaskForm
from wtforms import StringField, FloatField, SubmitField
from wtforms.validators import DataRequired, Length, NumberRange

class MoleculeForm(FlaskForm):
    name = StringField("Name", validators=[DataRequired(), Length(max=100)])
    formula = StringField("Formula", validators=[DataRequired(), Length(max=50)])
    molecular_weight = FloatField("Molecular Weight", validators=[DataRequired(), NumberRange(min=0)])
    submit = SubmitField("Add Molecule")
