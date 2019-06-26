"""
Creates a python file that can be run in Blender to generate a
representation of an object created by the CGinJulia package

'''Bl_Export(filename::String, points::AbstractVector{Point3H}, faces::AbstractVector{Real})'''
"""
function Bl_Export(filename::String,points,faces)

io = open(filename, "w")

Code::String = """bl_info = {
    "name": "Convex Hull Generation in Julia",
    "author": "Jackson Meade (NCSSM) and Jack Snoeyink (UNC Chapel Hill)",
    "version": (1, 0),
    "blender": (2, 80, 0),
    "location": "View3D > Add > Mesh > CGinJulia Convex Hull",
    "description": "Adds a new Mesh Object",
    "warning": "Please make sure you are familiar with the CGinJulia package before using this EXPERIMENTAL add-on",
    "wiki_url": "github.com/snoeyink/CGinJulia",
    "category": "Add Convex Hull",
}


import bpy
from bpy.types import Operator
from bpy.props import FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add
from mathutils import Vector


def add_object(self, context):
    scale_x = self.scale.x
    scale_y = self.scale.y
    scale_z = self.scale.z

    verts = ["""

for i = 1:length(points)
    Code = string(Code,string("""Vector((""",points[i][2],""" * scale_x""",""",""",points[i][3],""" * scale_y""",""",""",points[i][4],""" * scale_z""",""")),""","""
    """))
end

Code = string(Code,"""]
""","""    edges = []
""","""    faces = [""")

for a in faces
    Code = string(Code,string("""[""",string(a[1]),""",""",string(a[2]),""",""",string(a[3]),"""],"""))
end

Code = string(Code,"""]
""","""    mesh = bpy.data.meshes.new(name="Convex Hull of """,filename,"""")
    mesh.from_pydata(verts, edges, faces)
    # useful for development when the mesh may be invalid.
    # mesh.validate(verbose=True)
    object_data_add(context, mesh, operator=self)

class OBJECT_OT_add_object(Operator, AddObjectHelper):
    """,'\u0022','\u0022','\u0022',"""Generate the Mesh Object from Points and CGinJulia-Constructed Faces""",'\u0022','\u0022','\u0022',"""

    bl_idname = "mesh.add_object"
    bl_label = "Show Representation"
    bl_options = {'REGISTER', 'UNDO'}

    scale: FloatVectorProperty(
        name="scale",
        default=(1.0, 1.0, 1.0),
        subtype='TRANSLATION',
        description="scaling",
    )

    def execute(self, context):

        add_object(self, context)

        return {'FINISHED'}


# Registration

def add_object_button(self, context):
    self.layout.operator(
        OBJECT_OT_add_object.bl_idname,
        text="CGinJulia Convex Hull (""",filename,""")",
        icon='SURFACE_NSPHERE')

# This allows you to right click on a button and link to the manual
def add_object_manual_map():
    url_manual_prefix = "https://docs.blender.org/manual/en/dev/"
    url_manual_mapping = (
        ("bpy.ops.mesh.add_object", "editors/3dview/object"),
    )
    return url_manual_prefix, url_manual_mapping


def register():
    bpy.utils.register_class(OBJECT_OT_add_object)
    bpy.utils.register_manual_map(add_object_manual_map)
    bpy.types.VIEW3D_MT_mesh_add.append(add_object_button)


def unregister():
    bpy.utils.unregister_class(OBJECT_OT_add_object)
    bpy.utils.unregister_manual_map(add_object_manual_map)
    bpy.types.VIEW3D_MT_mesh_add.remove(add_object_button)


if __name__ == "__main__":
    register()
""")

write(io, Code)

close(io)

end
